#!/usr/bin/env Rscript
#
# write_inrich_intervals.R
#
# Per-trait step: reads locus files staged in the working directory for one
# trait and writes intervals_<trait_name>_for_INRICH.txt.
#
# SuSiE mode (--mode susie):
#   Reads trait<N>_locus_*.txt files. Interval = min/max bp of all SNPs in the
#   SuSiE window (already defined by LOD-drop in susie_finemap.R).
#
# Clump mode (--mode clump):
#   Reads trait<N>.clumped and trait<N>_loco.assoc.txt.
#   With --loddrop > 0: LOD-drop interval — walk left/right from each lead SNP
#   along the -log10(p) profile within the full LD neighbourhood (r² ≥ clump_r2,
#   within clump_kb), capped at maxdist. Mirrors ext_peak_sing() from mousegwas.
#   With --loddrop 0: fixed window of ±maxdist/2 around the lead SNP.
#
# Usage (called from WRITE_INRICH_INTERVALS Nextflow process):
#   Rscript write_inrich_intervals.R \
#       --trait_idx  3               \
#       --pheno_order phenotypes_order.txt

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--trait_idx",    type = "integer", required = TRUE,
    help = "1-based trait index")
parser$add_argument("--pheno_order",  required = TRUE,
    help = "phenotypes_order.txt from PREPARE_INPUT")
parser$add_argument("--mode",         default = "susie", choices = c("susie", "clump"),
    help = "Locus source: 'susie' or 'clump'")
parser$add_argument("--loddrop",      type = "double", default = 1.5,
    help = "LOD-drop threshold for interval expansion in clump mode (0 = fixed ±maxdist/2 window)")
parser$add_argument("--maxdist",      type = "double", default = 1e7,
    help = "Maximum distance in bp from peak for LOD-drop expansion (default: 10 Mb)")
parser$add_argument("--pval_type",    default = "p_wald",
    help = "p-value column in GEMMA assoc file (default: p_wald)")
parser$add_argument("--plink_prefix", default = NULL,
    help = "PLINK bfile prefix (bed/bim/fam without extension) for LD buddy calculation")
parser$add_argument("--clump_kb",     type = "double", default = 10000,
    help = "LD window in kb (must match the PLINK clump --clump-kb value)")
parser$add_argument("--clump_r2",     type = "double", default = 0.4,
    help = "r² threshold for LD buddy inclusion (must match --clump-r2)")
args <- parser$parse_args()

# ── Phenotype name lookup ──────────────────────────────────────────────────────
lines     <- readLines(args$pheno_order)
names_raw <- lines[-1]
names_raw <- trimws(gsub('^"|"$', "", names_raw))
names_raw <- names_raw[nzchar(names_raw)]

if (args$trait_idx > length(names_raw))
    stop(sprintf("trait_idx %d exceeds number of traits (%d)",
                 args$trait_idx, length(names_raw)))

trait_name <- gsub("[^A-Za-z0-9]", "_", names_raw[args$trait_idx])

# ── Helpers ────────────────────────────────────────────────────────────────────

# Returns all SNP IDs with r² ≥ clump_r2 to lead_snp within clump_kb, using
# PLINK --r2 --ld-snp. Includes the lead SNP itself. Returns character(0) on
# failure so the caller can fall back gracefully.
get_ld_buddies <- function(lead_snp, plink_prefix, clump_kb, clump_r2) {
    if (is.null(plink_prefix)) return(character(0))
    out_prefix <- tempfile(pattern = "ld_")
    cmd <- sprintf(
        "plink --bfile %s --r2 --ld-snp %s --ld-window 999999 --ld-window-kb %g --ld-window-r2 %g --out %s --silent",
        plink_prefix, lead_snp, clump_kb, clump_r2, out_prefix
    )
    ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    ld_file <- paste0(out_prefix, ".ld")
    if (ret != 0 || !file.exists(ld_file) || file.info(ld_file)$size == 0)
        return(character(0))
    ld <- tryCatch(fread(ld_file), error = function(e) data.table())
    unlink(c(ld_file, paste0(out_prefix, ".log"), paste0(out_prefix, ".nosex")),
           silent = TRUE)
    if (nrow(ld) == 0 || !"SNP_B" %in% names(ld)) return(character(0))
    unique(ld$SNP_B)
}

# LOD-drop interval — LD-neighbourhood-restricted, mirrors ext_peak_sing() in
# mousegwas. The walk is restricted to the lead SNP + all LD buddies (r² ≥
# clump_r2 within clump_kb), so it only crosses boundaries set by SNPs in the
# same LD block, not by unrelated low-p SNPs elsewhere on the chromosome.
#
# buddy_snps: character vector from get_ld_buddies() (may be empty).
# Returns a data.table with columns chr, minps, maxps.
lod_drop_interval <- function(assoc, lead_snp, buddy_snps, pval_col, loddrop, maxdist) {
    locus_snps  <- unique(c(lead_snp, buddy_snps))
    locus_assoc <- assoc[rs %in% locus_snps]

    peak_row <- locus_assoc[rs == lead_snp]
    if (nrow(peak_row) == 0) return(NULL)   # lead absent from assoc

    lead_chr  <- as.character(peak_row$chr[1])
    lead_bp   <- peak_row$ps[1]
    peak_pval <- peak_row[[pval_col]][1]

    if (is.na(peak_pval) || peak_pval <= 0 || !is.finite(-log10(peak_pval))) {
        return(data.table(chr   = lead_chr,
                          minps = pmax(1L, as.integer(lead_bp - maxdist / 2)),
                          maxps = as.integer(lead_bp + maxdist / 2)))
    }

    peak_logp <- -log10(peak_pval)
    lodstop   <- peak_logp - loddrop
    logp_vec  <- -log10(locus_assoc[[pval_col]])

    # Left boundary: rightmost locus SNP left of peak that drops below lodstop
    left_bp <- locus_assoc$ps[locus_assoc$ps < lead_bp & logp_vec < lodstop]
    minps   <- if (length(left_bp) > 0) max(lead_bp - maxdist, max(left_bp))
               else lead_bp - maxdist

    # Right boundary: leftmost locus SNP right of peak that drops below lodstop
    right_bp <- locus_assoc$ps[locus_assoc$ps > lead_bp & logp_vec < lodstop]
    maxps    <- if (length(right_bp) > 0) min(lead_bp + maxdist, min(right_bp))
                else lead_bp + maxdist

    data.table(chr   = lead_chr,
               minps = pmax(1L, as.integer(minps)),
               maxps = as.integer(maxps))
}

# ── Interval derivation ────────────────────────────────────────────────────────

if (args$mode == "susie") {
    # SuSiE mode: discover staged credible-set locus files for this trait
    locus_files <- Sys.glob(sprintf("trait%d_locus_*.txt", args$trait_idx))

    if (length(locus_files) == 0) {
        cat(sprintf("Trait %d (%s): no SuSiE locus files — skipping.\n",
                    args$trait_idx, trait_name))
        quit(status = 0)
    }

    loci <- rbindlist(lapply(locus_files, fread), fill = TRUE)

    if (nrow(loci) == 0) {
        cat(sprintf("Trait %d (%s): all locus files empty — skipping.\n",
                    args$trait_idx, trait_name))
        quit(status = 0)
    }

    # One interval per lead SNP: span of all SNPs in the SuSiE locus window
    intervals <- loci[, .(chr = chr[1], minps = min(bp), maxps = max(bp)),
                      by = lead_snp]

} else {
    # ── Clump mode ─────────────────────────────────────────────────────────────
    clump_file <- sprintf("trait%d.clumped", args$trait_idx)

    if (!file.exists(clump_file) || file.info(clump_file)$size == 0L) {
        cat(sprintf("Trait %d (%s): no clump file or no significant loci — skipping.\n",
                    args$trait_idx, trait_name))
        quit(status = 0)
    }

    clumped <- fread(clump_file, fill = TRUE)

    if (nrow(clumped) == 0 || !"SNP" %in% names(clumped)) {
        cat(sprintf("Trait %d (%s): clump file empty — skipping.\n",
                    args$trait_idx, trait_name))
        quit(status = 0)
    }

    lead_snps <- clumped[!is.na(SNP)]

    if (args$loddrop > 0) {
        # LOD-drop expansion using GEMMA assoc file + PLINK LD buddies
        assoc_file <- Sys.glob(sprintf("trait%d_loco.assoc.txt", args$trait_idx))

        if (length(assoc_file) == 0) {
            warning(sprintf("Trait %d: assoc file not found — falling back to fixed ±%.0f bp window.",
                            args$trait_idx, args$maxdist / 2))
            args$loddrop <- 0   # trigger fixed-window fallback below
        } else {
            assoc <- fread(assoc_file[1])

            intervals <- rbindlist(lapply(seq_len(nrow(lead_snps)), function(i) {
                row         <- lead_snps[i]
                buddy_snps  <- get_ld_buddies(row$SNP, args$plink_prefix,
                                              args$clump_kb, args$clump_r2)
                iv          <- lod_drop_interval(assoc, row$SNP, buddy_snps,
                                                 args$pval_type, args$loddrop,
                                                 args$maxdist)
                if (is.null(iv)) {
                    # Lead SNP filtered from assoc — use clump position as fallback
                    return(data.table(
                        lead_snp = row$SNP,
                        chr      = as.character(row$CHR),
                        minps    = pmax(1L, as.integer(row$BP - args$maxdist / 2)),
                        maxps    = as.integer(row$BP + args$maxdist / 2)
                    ))
                }
                cbind(data.table(lead_snp = row$SNP), iv)
            }), fill = TRUE)
        }
    }

    if (args$loddrop == 0) {
        # Fixed window: ±maxdist/2 around each lead SNP (mirrors loddrop=0 in mousegwas)
        half_win  <- as.integer(args$maxdist / 2)
        intervals <- lead_snps[, .(
            lead_snp = SNP,
            chr      = as.character(CHR[1]),
            minps    = pmax(1L, as.integer(BP[1] - half_win)),
            maxps    = as.integer(BP[1] + half_win)
        ), by = SNP][, SNP := NULL]
    }
}

out_file <- paste0("intervals_", trait_name, "_for_INRICH.txt")
write.table(intervals[, .(chr, minps, maxps)],
            out_file, sep = "\t", col.names = FALSE,
            row.names = FALSE, quote = FALSE)

cat(sprintf("Trait %d (%s): %d interval(s) written to %s\n",
            args$trait_idx, trait_name, nrow(intervals), out_file))
