#!/usr/bin/env Rscript
#
# susie_finemap.R
#
# Runs SuSiE fine-mapping for all loci identified by PLINK clumping for a
# single trait. Loci are defined by the lead SNPs in the clump file, expanded
# by --window bp on each side. Uses individual-level genotypes via BEDMatrix.
#
# Outputs (written to --outdir):
#   trait{N}_locus_{lead_snp}.txt  — per-SNP PIP and credible set membership

library(argparse)
library(data.table)
library(BEDMatrix)
library(susieR)

parser <- ArgumentParser()
parser$add_argument("--bfile",     required = TRUE,  help = "PLINK BED prefix")
parser$add_argument("--pheno",     required = TRUE,  help = "GEMMA-format phenotype file (no header, no IID)")
parser$add_argument("--trait",     required = TRUE,  type = "integer", help = "Phenotype column index")
parser$add_argument("--gwas",      required = TRUE,  help = "GEMMA LOCO .assoc.txt file for this trait")
parser$add_argument("--clump",     required = TRUE,  help = "PLINK .clumped file for this trait")
parser$add_argument("--threshold", required = TRUE,  help = "File containing significance threshold (single value)")
parser$add_argument("--window",    default  = 1e6,   type = "double",
                    help = "Max half-window size in bp (cap for LOD-drop, or fixed window when --loddrop 0; default 1Mb)")
parser$add_argument("--loddrop",   default  = 1.5,   type = "double",
                    help = "LOD-drop threshold for interval expansion: extend from peak until -log10(p) drops by this much (0 = use fixed --window; default 1.5)")
parser$add_argument("--pval_type", default  = "p_wald",
                    help = "P-value column to use for LOD-drop (default: p_wald)")
parser$add_argument("--L",         default  = 10,    type = "integer",
                    help = "Max number of causal signals per locus (default 10)")
parser$add_argument("--outdir",    default  = ".")
args <- parser$parse_args()

# ── Load threshold ────────────────────────────────────────────────────────────
pthresh <- as.numeric(readLines(args$threshold, n = 1L))
cat(sprintf("Significance threshold: %.2e\n", pthresh))

# ── Load clump file ───────────────────────────────────────────────────────────
clump <- fread(args$clump)
if (nrow(clump) == 0) {
    cat(sprintf("Trait %d: no significant loci in clump file, nothing to fine-map.\n", args$trait))
    quit(status = 0)
}

cat(sprintf("Trait %d: %d locus/loci from clump file.\n", args$trait, nrow(clump)))

# ── Load BED/BIM/FAM ──────────────────────────────────────────────────────────
bim <- fread(paste0(args$bfile, ".bim"), header = FALSE,
             col.names = c("chr", "snp", "cm", "bp", "a1", "a2"))
fam <- fread(paste0(args$bfile, ".fam"), header = FALSE)
n   <- nrow(fam)
p   <- nrow(bim)
bed <- BEDMatrix(paste0(args$bfile, ".bed"), n = n, p = p)

# ── Load phenotype ────────────────────────────────────────────────────────────
pheno_mat <- as.matrix(fread(args$pheno, header = FALSE,
                             na.strings = c("NA", "N/A", "na", "")))
y <- as.numeric(pheno_mat[, args$trait])

# ── Load GWAS results ─────────────────────────────────────────────────────────
gwas <- fread(args$gwas)
if (!args$pval_type %in% names(gwas))
    stop(sprintf("Column '%s' not found in GWAS file.", args$pval_type))
gwas[, logp := -log10(get(args$pval_type))]

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

# ── Run SuSiE per locus ───────────────────────────────────────────────────────
for (i in seq_len(nrow(clump))) {
    lead_snp <- clump$SNP[i]
    lead_chr <- clump$CHR[i]
    lead_bp  <- clump$BP[i]

    if (args$loddrop > 0) {
        gwas_chr  <- gwas[as.character(chr) == as.character(lead_chr)]
        peak_logp <- gwas_chr[ps == lead_bp, logp]
        if (length(peak_logp) == 0) peak_logp <- -log10(clump$P[i])
        lodstop   <- peak_logp - args$loddrop
        win_start <- max(1L, max(c(lead_bp - args$window,
                                   gwas_chr[ps < lead_bp & logp < lodstop, ps])))
        win_end   <-         min(c(lead_bp + args$window,
                                   gwas_chr[ps > lead_bp & logp < lodstop, ps]))
    } else {
        # loddrop=0: span of all clumped SNPs capped at ±window (matches mousegwas behavior)
        sp2_raw  <- clump$SP2[i]
        sp2_snps <- if (is.na(sp2_raw) || trimws(sp2_raw) %in% c("NONE", "")) {
            character(0)
        } else {
            trimws(sub("\\([^)]+\\)", "", strsplit(sp2_raw, ",")[[1]]))
        }
        clump_bps <- bim$bp[bim$snp %in% c(lead_snp, sp2_snps)]
        if (length(clump_bps) == 0) clump_bps <- lead_bp
        win_start <- max(1L, max(min(clump_bps), lead_bp - args$window))
        win_end   <-         min(max(clump_bps), lead_bp + args$window)
    }

    snp_idx <- which(bim$chr == lead_chr &
                     bim$bp  >= win_start &
                     bim$bp  <= win_end)

    if (length(snp_idx) < 2) {
        cat(sprintf("  Locus %s: fewer than 2 SNPs in window, skipping.\n", lead_snp))
        next
    }

    interval_type <- if (args$loddrop > 0) sprintf("LOD-drop %.1f", args$loddrop) else "fixed window"
    cat(sprintf("  Locus %s | chr%s:%d-%d | %d SNPs [%s]\n",
                lead_snp, lead_chr, win_start, win_end, length(snp_idx), interval_type))

    # Load, impute missing with column mean, and standardise
    X <- as.matrix(bed[, snp_idx])
    colnames(X) <- bim$snp[snp_idx]

    col_means <- colMeans(X, na.rm = TRUE)
    for (j in seq_len(ncol(X))) {
        miss <- is.na(X[, j])
        if (any(miss)) X[miss, j] <- col_means[j]
    }

    X    <- scale(X)
    keep <- !is.nan(X[1, ])   # drop monomorphic SNPs
    X       <- X[, keep, drop = FALSE]
    snp_idx <- snp_idx[keep]

    if (ncol(X) < 2) {
        cat(sprintf("  Locus %s: fewer than 2 polymorphic SNPs after QC, skipping.\n", lead_snp))
        next
    }

    obs   <- !is.na(y)
    X_fit <- X[obs, , drop = FALSE]
    y_fit <- y[obs]

    fit <- susie(X_fit, y_fit,
                 L                         = args$L,
                 estimate_residual_variance = TRUE,
                 estimate_prior_variance    = TRUE,
                 verbose                    = FALSE)

    results <- data.table(
        trait    = args$trait,
        lead_snp = lead_snp,
        snp      = bim$snp[snp_idx][keep],
        chr      = bim$chr[snp_idx][keep],
        bp       = bim$bp[snp_idx][keep],
        pip      = fit$pip,
        cs       = NA_integer_
    )

    n_cs <- 0
    if (!is.null(fit$sets$cs)) {
        n_cs <- length(fit$sets$cs)
        for (cs_i in seq_len(n_cs))
            results[fit$sets$cs[[cs_i]], cs := cs_i]
    }

    out_file <- file.path(args$outdir,
                          sprintf("trait%d_locus_%s.txt", args$trait, lead_snp))
    fwrite(results, out_file, sep = "\t")
    cat(sprintf("    -> %d credible set(s). Saved to %s\n", n_cs, out_file))
}
