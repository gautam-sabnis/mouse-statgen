#!/usr/bin/env Rscript
#
# plot_inrich.R
#
# Collects INRICH results across all labels (per-trait, per-group, per-cluster)
# and gene-set types, then produces one heatmap PDF per gene-set type showing
# adjusted p-values (rows = labels, columns = enriched terms).
#
# Adapted from mousegwas/exec/plot_INRICH_pvalues.R.
#
# Usage:
#   Rscript plot_inrich.R          \
#       --files  *.out.inrich      \
#       --yaml   run.yaml          \
#       --name   <run_name>        \
#       --outdir .                 \
#       --minpval 0.05

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ComplexHeatmap)
    library(yaml)
    library(RColorBrewer)
    library(argparse)
    library(grid)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

parser <- ArgumentParser()
parser$add_argument("--files",   "-f", nargs = "+", required = TRUE,
    help = "List of INRICH .out.inrich result files")
parser$add_argument("--yaml",    "-y", required = TRUE,
    help = "Run YAML config (used to map trait names to paper names and groups)")
parser$add_argument("--name",    required = TRUE, help = "Run name")
parser$add_argument("--outdir",  "-o", default = ".",
    help = "Output directory for heatmap PDFs")
parser$add_argument("--minpval", "-p", type = "double", default = 0.05,
    help = "Only show gene-set terms with adjusted p-value <= this threshold in at least one label")
args <- parser$parse_args()

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# ── Parse one INRICH output file ─────────────────────────────────────────────
# INRICH writes lines containing "_O1" for the observed statistics.
parse_inrich <- function(fname) {
    con <- pipe(paste0('grep "_O1" ', shQuote(fname)))
    tbl <- tryCatch(
        read.table(con, sep = "\t", header = TRUE, quote = ""),
        error = function(e) NULL
    )
    try(close(con), silent = TRUE)
    tbl
}

# ── Extract label and gene-set from filename ──────────────────────────────────
# Expected basename: intervals_<label>_groups_<geneset>.out.inrich
extract_meta <- function(fname) {
    bn    <- sub("\\.out\\.inrich$", "", basename(fname))
    label <- sub("^intervals_(.*)_groups_.*$", "\\1", bn)
    gset  <- sub("^intervals_.*_groups_(.*)$", "\\1", bn)
    list(label = label, geneset = gset)
}

# ── Load YAML for trait name → paper name + group mapping ────────────────────
yamin <- yaml.load_file(args$yaml)

# Build lookup: internal name → (PaperName, Group)
pnames <- tibble(
    internal = character(0),
    PaperName = character(0),
    Group = character(0)
)
for (n in names(yamin$phenotypes)) {
    pnames <- bind_rows(pnames, tibble(
        internal  = n,
        PaperName = yamin$phenotypes[[n]]$papername %||% n,
        Group     = yamin$phenotypes[[n]]$group %||% "NoGroup"
    ))
}

# Also build a lookup keyed by sanitised internal name (matching our label format)
pnames$label_key <- gsub("[^A-Za-z0-9]", "_", pnames$internal)

# Group colour palette
groupsOrder <- if (length(yamin$groups)) c(yamin$groups, "General") else
               c(unique(pnames$Group), "General")
grpcol <- rep(brewer.pal(8, "Accent"), ceiling(length(groupsOrder) / 8))[seq_along(groupsOrder)]
names(grpcol) <- groupsOrder

# ── Collect all results ───────────────────────────────────────────────────────
all_results <- NULL
for (f in args$files) {
    meta <- extract_meta(f)
    tbl  <- parse_inrich(f)
    if (is.null(tbl) || nrow(tbl) == 0) next
    tbl$label   <- meta$label
    tbl$geneset <- meta$geneset
    all_results <- bind_rows(all_results, tbl)
}

if (is.null(all_results) || nrow(all_results) == 0) {
    cat("No INRICH results to plot.\n")
    quit(status = 0)
}

all_results <- as_tibble(all_results)
all_results <- all_results %>%
    mutate(logpval = -log10(PCORR))

# ── Resolve display names for labels ─────────────────────────────────────────
# Per-trait labels map to YAML paper names; group/cluster labels are used as-is.
resolve_display <- function(label) {
    # Try matching to YAML internal name
    idx <- match(label, pnames$label_key)
    if (!is.na(idx)) return(pnames$PaperName[idx])
    # Try matching the raw internal name
    idx2 <- match(label, gsub("[^A-Za-z0-9]", "_", pnames$PaperName))
    if (!is.na(idx2)) return(pnames$PaperName[idx2])
    # Group / cluster labels: pretty-print
    gsub("_", " ", label)
}

resolve_group <- function(label) {
    idx <- match(label, pnames$label_key)
    if (!is.na(idx)) return(pnames$Group[idx])
    if (grepl("^cluster_", label)) return("Cluster")
    "General"
}

unique_labels <- unique(all_results$label)
label_meta <- tibble(
    label       = unique_labels,
    display     = sapply(unique_labels, resolve_display),
    group       = sapply(unique_labels, resolve_group)
) %>%
    mutate(color = grpcol[group] %||% grpcol["General"])

# ── Plot one heatmap per gene-set type ────────────────────────────────────────
plot_heatmap <- function(tbl, title) {
    # Filter to terms significant in at least one label
    sig_terms <- unique(tbl$TARGET[tbl$PCORR <= args$minpval])
    if (length(sig_terms) == 0) {
        cat(sprintf("  %s: no significant terms at p <= %.3f\n", title, args$minpval))
        return(invisible(NULL))
    }
    tbl <- tbl %>% filter(TARGET %in% sig_terms)

    mat <- tbl %>%
        dplyr::select(label, TARGET, logpval) %>%
        pivot_wider(names_from = TARGET, values_from = logpval,
                    values_fn = max, values_fill = 0) %>%
        column_to_rownames("label") %>%
        as.matrix()

    # Align row metadata
    rm <- label_meta[match(rownames(mat), label_meta$label), ]
    row_ann <- rowAnnotation(
        Group = rm$group,
        col   = list(Group = grpcol[unique(rm$group)]),
        show_legend = FALSE
    )
    lgp_at  <- seq(0, ceiling(max(mat, na.rm = TRUE)))
    lgp_lbl <- formatC(10^(-lgp_at), format = "e", digits = 0)

    h <- Heatmap(
        mat,
        name              = "adj. p-value",
        row_labels        = rm$display,
        column_title      = title,
        row_title         = "Labels",
        heatmap_legend_param = list(at = lgp_at, labels = lgp_lbl),
        right_annotation  = row_ann,
        column_names_gp   = grid::gpar(fontsize = 7),
        row_names_gp      = grid::gpar(fontsize = 8)
    )
    draw(h, padding = unit(c(1.5, 0.1, 0.1, 0.2), "inches"))
}

for (gset in unique(all_results$geneset)) {
    tbl_g  <- filter(all_results, geneset == gset)
    outpdf <- file.path(args$outdir,
                        paste0("inrich_heatmap_", args$name, "_", gset, ".pdf"))
    cairo_pdf(outpdf, width = 8.25, height = 11.25)
    plot_heatmap(tbl_g, gset)
    dev.off()
    cat(sprintf("Written: %s\n", outpdf))
}
