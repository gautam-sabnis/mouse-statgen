#!/usr/bin/env Rscript
#
# write_inrich_intervals_combined.R
#
# Runs once after CLUSTER_PEAKS. Aggregates locus files across all traits to
# build interval files at three levels:
#   - All phenotypes combined
#   - Per phenotype group   (from trait_groups TSV)
#   - Per k-means cluster   (from cluster_assignments_*.csv)
#
# SuSiE mode (--mode susie):
#   Reads trait*_locus_*.txt. Interval = min/max bp across all SNPs in the window.
#
# Clump mode (--mode clump):
#   Reads trait*.clumped and trait*_loco.assoc.txt files.
#   With --loddrop > 0: LOD-drop intervals within the full LD neighbourhood
#   (r² ≥ clump_r2 within clump_kb), capped at maxdist.
#   With --loddrop 0: fixed ±maxdist/2 window around each lead SNP.
#   Intervals are unioned per lead_snp across all traits that share that locus.

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--pheno_order",  required = TRUE)
parser$add_argument("--trait_groups", required = TRUE,
    help = "trait_groups_<name>.tsv: two-column TSV (trait_name, group), with header")
parser$add_argument("--name",         required = TRUE, help = "Run name (for messages only)")
parser$add_argument("--mode",         default = "susie", choices = c("susie", "clump"),
    help = "Locus source: 'susie' (credible set files) or 'clump' (PLINK .clumped files)")
parser$add_argument("--loddrop",      type = "double", default = 1.5,
    help = "LOD-drop threshold for interval expansion in clump mode (0 = fixed window)")
parser$add_argument("--maxdist",      type = "double", default = 1e7,
    help = "Max distance in bp from peak for LOD-drop expansion (default: 10 Mb)")
parser$add_argument("--pval_type",    default = "p_wald",
    help = "p-value column in GEMMA assoc files (default: p_wald)")
parser$add_argument("--plink_prefix", default = NULL,
    help = "PLINK bfile prefix for LD buddy calculation")
parser$add_argument("--clump_kb",     type = "double", default = 10000,
    help = "LD window in kb (must match --clump-kb)")
parser$add_argument("--clump_r2",     type = "double", default = 0.4,
    help = "r² threshold for LD buddy inclusion (must match --clump-r2)")
args <- parser$parse_args()

# ── Write helper: one row per unique lead_snp ─────────────────────────────────
# dt must have columns: lead_snp, chr, bp_min, bp_max
# For each lead_snp, take min(bp_min) and max(bp_max) to union across traits.
write_intervals <- function(dt, label) {
    if (nrow(dt) == 0) {
        cat(sprintf("  %s: no loci, skipping.\n", label))
        return(invisible(NULL))
    }
    intervals <- dt[, .(chr = chr[1], minps = min(bp_min), maxps = max(bp_max)),
                    by = lead_snp]
    fname <- paste0("intervals_", label, "_for_INRICH.txt")
    write.table(intervals[, .(chr, minps, maxps)],
                fname, sep = "\t", col.names = FALSE,
                row.names = FALSE, quote = FALSE)
    cat(sprintf("  %s: %d intervals written.\n", label, nrow(intervals)))
}

# ── Helpers ────────────────────────────────────────────────────────────────────

# Returns all SNP IDs with r² ≥ clump_r2 to lead_snp within clump_kb using
# PLINK --r2 --ld-snp. Results are cached to avoid redundant PLINK calls when
# the same lead SNP appears across multiple traits.
ld_cache <- list()

get_ld_buddies <- function(lead_snp, plink_prefix, clump_kb, clump_r2) {
    if (is.null(plink_prefix)) return(character(0))
    if (!is.null(ld_cache[[lead_snp]])) return(ld_cache[[lead_snp]])
    out_prefix <- tempfile(pattern = "ld_")
    cmd <- sprintf(
        "plink --bfile %s --r2 --ld-snp %s --ld-window 999999 --ld-window-kb %g --ld-window-r2 %g --out %s --silent",
        plink_prefix, lead_snp, clump_kb, clump_r2, out_prefix
    )
    ret <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    ld_file <- paste0(out_prefix, ".ld")
    result <- character(0)
    if (ret == 0 && file.exists(ld_file) && file.info(ld_file)$size > 0) {
        ld <- tryCatch(fread(ld_file), error = function(e) data.table())
        if (nrow(ld) > 0 && "SNP_B" %in% names(ld))
            result <- unique(ld$SNP_B)
    }
    unlink(c(ld_file, paste0(out_prefix, ".log"), paste0(out_prefix, ".nosex")),
           silent = TRUE)
    ld_cache[[lead_snp]] <<- result
    result
}

# LOD-drop interval — LD-neighbourhood-restricted, mirrors ext_peak_sing() in
# mousegwas. Walk is restricted to the lead SNP + all r²-LD buddies so it only
# crosses boundaries set by SNPs in the same LD block.
lod_drop_interval <- function(assoc, lead_snp, buddy_snps, pval_col, loddrop, maxdist) {
    locus_snps  <- unique(c(lead_snp, buddy_snps))
    locus_assoc <- assoc[rs %in% locus_snps]

    peak_row <- locus_assoc[rs == lead_snp]
    if (nrow(peak_row) == 0) return(NULL)

    lead_chr  <- as.character(peak_row$chr[1])
    lead_bp   <- peak_row$ps[1]
    peak_pval <- peak_row[[pval_col]][1]

    if (is.na(peak_pval) || peak_pval <= 0 || !is.finite(-log10(peak_pval))) {
        return(list(chr    = lead_chr,
                    bp_min = pmax(1L, as.integer(lead_bp - maxdist / 2)),
                    bp_max = as.integer(lead_bp + maxdist / 2)))
    }

    peak_logp <- -log10(peak_pval)
    lodstop   <- peak_logp - loddrop
    logp_vec  <- -log10(locus_assoc[[pval_col]])

    left_bp <- locus_assoc$ps[locus_assoc$ps < lead_bp & logp_vec < lodstop]
    minps   <- if (length(left_bp) > 0) max(lead_bp - maxdist, max(left_bp))
               else lead_bp - maxdist

    right_bp <- locus_assoc$ps[locus_assoc$ps > lead_bp & logp_vec < lodstop]
    maxps    <- if (length(right_bp) > 0) min(lead_bp + maxdist, min(right_bp))
                else lead_bp + maxdist

    list(chr    = lead_chr,
         bp_min = pmax(1L, as.integer(minps)),
         bp_max = as.integer(maxps))
}

# ── Load phenotype names ──────────────────────────────────────────────────────
lines     <- readLines(args$pheno_order)
names_raw <- trimws(gsub('^"|"$', "", lines[-1]))
names_raw <- names_raw[nzchar(names_raw)]

# ── Load locus data ───────────────────────────────────────────────────────────
if (args$mode == "susie") {
    locus_files <- Sys.glob("trait*_locus_*.txt")

    if (length(locus_files) == 0) {
        cat("No SuSiE locus files found — no combined interval files written.\n")
        quit(status = 0)
    }

    all_loci <- rbindlist(lapply(locus_files, function(f) {
        dt <- fread(f, fill = TRUE)
        dt[, trait_idx := as.integer(sub("trait(\\d+)_locus_.*", "\\1", basename(f)))]
        dt
    }), fill = TRUE)

    all_loci[, trait_name := gsub("[^A-Za-z0-9]", "_", names_raw[trait_idx])]
    # SuSiE: interval = span of all SNPs in the window
    all_loci[, `:=`(bp_min = bp, bp_max = bp)]

} else {
    # ── Clump mode ─────────────────────────────────────────────────────────────
    clump_files <- Sys.glob("trait*.clumped")

    if (length(clump_files) == 0) {
        cat("No clump files found — no combined interval files written.\n")
        quit(status = 0)
    }

    all_loci <- rbindlist(lapply(clump_files, function(f) {
        dt <- tryCatch(fread(f, fill = TRUE), error = function(e) data.table())
        if (nrow(dt) == 0 || !"SNP" %in% names(dt)) return(data.table())
        idx <- as.integer(sub("trait(\\d+)\\.clumped", "\\1", basename(f)))
        dt[!is.na(SNP), .(lead_snp = SNP, trait_idx = idx,
                           chr = as.character(CHR), bp = BP)]
    }), fill = TRUE)

    if (nrow(all_loci) == 0) {
        cat("No significant clump loci found — no combined interval files written.\n")
        quit(status = 0)
    }

    all_loci[, trait_name := gsub("[^A-Za-z0-9]", "_", names_raw[trait_idx])]

    # ── Compute intervals via LOD-drop or fixed window ─────────────────────────
    if (args$loddrop > 0) {
        assoc_files <- Sys.glob("trait*_loco.assoc.txt")

        if (length(assoc_files) == 0) {
            warning("No assoc files found — falling back to fixed ±maxdist/2 window.")
            args$loddrop <- 0
        } else {
            # Build lookup: trait_idx -> assoc data.table
            assoc_by_trait <- setNames(
                lapply(assoc_files, fread),
                as.integer(sub("trait(\\d+)_loco\\.assoc\\.txt", "\\1",
                               basename(assoc_files)))
            )

            interval_rows <- rbindlist(lapply(seq_len(nrow(all_loci)), function(i) {
                row         <- all_loci[i]
                assoc       <- assoc_by_trait[[as.character(row$trait_idx)]]
                buddy_snps  <- get_ld_buddies(row$lead_snp, args$plink_prefix,
                                              args$clump_kb, args$clump_r2)
                fallback <- data.table(lead_snp  = row$lead_snp,
                                       trait_idx = row$trait_idx,
                                       chr       = row$chr,
                                       bp_min    = pmax(1L, as.integer(row$bp - args$maxdist / 2)),
                                       bp_max    = as.integer(row$bp + args$maxdist / 2))
                if (is.null(assoc)) return(fallback)
                iv <- lod_drop_interval(assoc, row$lead_snp, buddy_snps,
                                        args$pval_type, args$loddrop, args$maxdist)
                if (is.null(iv)) return(fallback)
                data.table(lead_snp  = row$lead_snp,
                           trait_idx = row$trait_idx,
                           chr       = iv$chr,
                           bp_min    = iv$bp_min,
                           bp_max    = iv$bp_max)
            }), fill = TRUE)

            all_loci <- merge(all_loci, interval_rows,
                              by = c("lead_snp", "trait_idx", "chr"), all.x = TRUE)
        }
    }

    if (args$loddrop == 0) {
        half_win <- as.integer(args$maxdist / 2)
        all_loci[, `:=`(bp_min = pmax(1L, as.integer(bp - half_win)),
                         bp_max = as.integer(bp + half_win))]
    }
}

cat(sprintf("Loaded %d loci across %d traits.\n",
    uniqueN(all_loci$lead_snp), uniqueN(all_loci$trait_idx)))

# ── 1. All phenotypes combined ────────────────────────────────────────────────
cat("Writing combined interval files...\n")
write_intervals(all_loci, "All_Phenotypes")

# ── 2. Per phenotype group ────────────────────────────────────────────────────
tgroups <- fread(args$trait_groups)
colnames(tgroups) <- c("trait_name", "group")

all_loci_grp <- merge(all_loci, tgroups, by = "trait_name", all.x = TRUE)

for (grp in sort(unique(tgroups$group))) {
    safe <- gsub("[^A-Za-z0-9]", "_", grp)
    write_intervals(all_loci_grp[group == grp], safe)
}

# ── 3. Per k-means cluster ────────────────────────────────────────────────────
cluster_files <- Sys.glob("cluster_assignments_*.csv")

if (length(cluster_files) == 0) {
    cat("No cluster_assignments file found — skipping per-cluster intervals.\n")
} else {
    cluster_assign <- fread(cluster_files[1])
    all_loci_cl <- merge(all_loci, cluster_assign[, .(lead_snp, cluster)],
                         by = "lead_snp", all.x = FALSE)
    for (k in sort(unique(all_loci_cl$cluster))) {
        write_intervals(all_loci_cl[cluster == k], paste0("cluster_", k))
    }
}

cat("Done.\n")
