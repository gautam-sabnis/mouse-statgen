#!/usr/bin/env Rscript
#
# plot_pve.R
#
# Plots heritability (h2 / PVE) estimates with standard errors per phenotype,
# colored by phenotype group.
#
# Accepts output from any of the three estimation methods in the pipeline.
# Column names are normalized internally so plot aesthetics are identical
# regardless of method:
#
#   GEMMA  heritability_gemma_<name>.csv  →  PVE, PVESE  (renamed to h2, se)
#   GCTA   heritability_gcta_<name>.csv   →  h2, se
#   LDAK   heritability_<name>.csv        →  h2, se
#
# Outputs:
#   PVE_plot_<method>_<name>.pdf

# ── Sourcing ───────────────────────────────────────────────────────────────────
# LOCAL <- TRUE  : run interactively (RStudio / Rscript from project root);
#                  argparse is skipped, paths taken from the args list below.
# LOCAL <- FALSE : run via Rscript on HPC / inside Nextflow; argparse is used.
LOCAL <- FALSE

if (LOCAL) {
    SRC <- "src"
} else {
    SRC <- dirname(normalizePath(
        Filter(\(x) grepl("--file=", x), commandArgs(trailingOnly = FALSE)) |>
            sub("--file=", "", x = _)
    ))
}
source(file.path(SRC, "aesthetics.R"))
source(file.path(SRC, "postprocess_setup.R"))

# ── Packages ───────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(argparse)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
# LOCAL mode: update paths below, then set LOCAL <- TRUE above.
if (LOCAL) {
    args <- list(
        input    = "output/heritability/heritability_gemma_anxfb_test.csv",
        method   = "gemma",
        name     = "anxfb_test",
        yaml     = "data/anxfb.yaml",
        outdir   = "output/plots",
        minherit = 0,
        palette  = "nejm"
    )
} else {
    parser <- ArgumentParser()
    parser$add_argument("--input",  "-i", required = TRUE,
        help = "Heritability CSV from COLLECT_HERITABILITY, COLLECT_HERITABILITY_GCTA, or HERITABILITY_LDAK")
    parser$add_argument("--method", "-m", required = TRUE,
        choices = c("gemma", "gcta", "ldak"),
        help = "Estimation method — determines column mapping and output filename")
    parser$add_argument("--name", "-n", required = TRUE,
        help = "Run name (params.name from Nextflow, e.g. anxfb_test) — used in output filename")
    parser$add_argument("--yaml", "-y", required = TRUE,
        help = "Pipeline YAML with phenotype groups and paper names")
    parser$add_argument("--outdir", "-o", default = ".",
        help = "Output directory (default: current dir)")
    parser$add_argument("--minherit", type = "double", default = 0,
        help = "h2 threshold: bars below this are shown at reduced opacity (default: 0)")
    parser$add_argument("--palette", default = "nejm",
        choices = c("nejm", "npg", "bmj", "jco", "lancet", "jama"),
        help = "ggsci color palette for phenotype groups (default: nejm)")
    args <- parser$parse_args()
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# ── Phenotype metadata ─────────────────────────────────────────────────────────
cfg              <- load_phenotype_config(args$yaml)
colors           <- assign_group_colors(cfg$pnames, cfg$groups_order, palette = args$palette)
pnames           <- colors$pnames
rownames(pnames) <- cfg$pheno_names
grpcol           <- colors$grpcol

# ── Load and normalize heritability estimates ──────────────────────────────────
h2_raw <- readr::read_csv(args$input, show_col_types = FALSE)

# GEMMA uses PVE / PVESE; GCTA and LDAK already use h2 / se
if (args$method == "gemma") {
    h2_raw <- dplyr::rename(h2_raw, h2 = PVE, se = PVESE)
}

# Join with phenotype metadata
h2 <- dplyr::left_join(
    h2_raw,
    tibble::as_tibble(pnames, rownames = "trait"),
    by = "trait"
)

# ── Plot ───────────────────────────────────────────────────────────────────────
pvh <- if (nrow(h2) > 40) FIG$panel_height * 2 else FIG$panel_height
out <- file.path(args$outdir, paste0("PVE_plot_", args$method, "_", args$name, ".pdf"))

p <- ggplot(h2, aes(x     = reorder(PaperName, -h2),
                    y     = h2,
                    fill  = Group,
                    alpha = h2 > args$minherit)) +
    geom_bar(color = "black", stat = "identity") +
    geom_errorbar(aes(ymin = h2 - se, ymax = h2 + se), width = 0.2) +
    scale_fill_manual(values = grpcol) +
    scale_alpha_manual(values = c(`FALSE` = 0.4, `TRUE` = 1), guide = "none") +
    xlab("Phenotype") +
    ylab(bquote(h^2 ~ "(±SE)")) +
    coord_flip() +
    theme_pve()

save_plot(p, path = out, width = FIG$half_width * 1.5, height = pvh)
message("Written: ", out)
