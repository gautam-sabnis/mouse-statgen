#!/usr/bin/env Rscript
#
# collect_genetic_corr.R
#
# Assembles pairwise GCTA --reml-bivar outputs into symmetric n x n matrices
# for genetic correlation (rG) and LRT p-value (Pval).
#
# Each input file is named rg_i_j.hsq (i > j, lower-triangular pairs).
#
# Outputs:
#   genetic_corr_<name>.csv   -- rG matrix  (diagonal = 1,  NA = failed pairs)
#   genetic_pval_<name>.csv   -- Pval matrix (diagonal = NA, NA = failed pairs)

library(argparse)
library(data.table)

parser <- ArgumentParser()
parser$add_argument("--hsq_dir",     required = TRUE,
    help = "Directory containing rg_i_j.hsq files")
parser$add_argument("--pheno_order", required = TRUE,
    help = "phenotypes_order.txt from PREPARE_INPUT (header + one trait per line)")
parser$add_argument("--name",        required = TRUE,
    help = "Run name (used in output filenames)")
parser$add_argument("--outdir",      default  = ".",
    help = "Output directory (default: current dir)")
args <- parser$parse_args()

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

# ── 1. Trait names (1-based index → name) ─────────────────────────────────────
trait_names <- read.csv(args$pheno_order, header = FALSE, skip = 1)$V1

# ── 2. Parse all .hsq files ───────────────────────────────────────────────────
files <- list.files(args$hsq_dir,
                    pattern    = "^rg_\\d+_\\d+\\.hsq$",
                    full.names = TRUE)
if (length(files) == 0)
    stop("No rg_i_j.hsq files found in: ", args$hsq_dir)

cat(sprintf("Parsing %d .hsq files...\n", length(files)))

pairs <- lapply(files, function(f) {
    bn  <- sub("\\.hsq$", "", basename(f))
    idx <- as.integer(strsplit(sub("^rg_", "", bn), "_")[[1]])
    i   <- idx[1]; j <- idx[2]

    lines <- readLines(f, warn = FALSE)

    rg_line   <- lines[grepl("^rG\t", lines)]
    pval_line <- lines[grepl("^Pval\t", lines)]

    if (length(rg_line) == 0) {
        warning(sprintf("No rG in %s — skipping (convergence failure?)", basename(f)))
        return(NULL)
    }

    rg_vals <- strsplit(trimws(rg_line[1]),   "\t")[[1]]
    rg      <- as.numeric(rg_vals[2])
    se      <- as.numeric(rg_vals[3])

    pval <- if (length(pval_line) > 0) {
        pv <- strsplit(trimws(pval_line[1]), "\t")[[1]][2]
        as.numeric(gsub("\\s*\\(.*\\)", "", pv))   # strip "(one-tailed test)" if present
    } else {
        NA_real_
    }

    list(i = i, j = j, rg = rg, se = se, pval = pval)
})
pairs <- Filter(Negate(is.null), pairs)
cat(sprintf("  %d pairs parsed successfully, %d failed.\n",
            length(pairs), length(files) - length(pairs)))

# ── 3. Build symmetric matrices ───────────────────────────────────────────────
all_idx  <- sort(unique(c(sapply(pairs, `[[`, "i"), sapply(pairs, `[[`, "j"))))
labels   <- trait_names[all_idx]
n        <- length(all_idx)
idx_map  <- setNames(seq_len(n), all_idx)

rg_mat   <- matrix(NA_real_, n, n, dimnames = list(labels, labels))
pval_mat <- matrix(NA_real_, n, n, dimnames = list(labels, labels))
diag(rg_mat) <- 1   # self-correlation = 1; p-value diagonal stays NA

for (p in pairs) {
    ri <- idx_map[as.character(p$i)]
    ci <- idx_map[as.character(p$j)]
    rg_mat[ri, ci]   <- p$rg;   rg_mat[ci, ri]   <- p$rg
    pval_mat[ri, ci] <- p$pval; pval_mat[ci, ri] <- p$pval
}

# ── 4. Write outputs ──────────────────────────────────────────────────────────
out_rg   <- file.path(args$outdir, sprintf("genetic_corr_%s.csv",  args$name))
out_pval <- file.path(args$outdir, sprintf("genetic_pval_%s.csv",  args$name))

fwrite(as.data.table(rg_mat,   keep.rownames = "trait"), out_rg)
fwrite(as.data.table(pval_mat, keep.rownames = "trait"), out_pval)

cat(sprintf("Saved: %s\n", out_rg))
cat(sprintf("Saved: %s\n", out_pval))
cat(sprintf("Matrix dimensions: %d x %d traits\n", n, n))
