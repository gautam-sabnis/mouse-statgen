#!/usr/bin/env Rscript
#
# plot_manhattan.R
#
# Plots Manhattan plots for individual traits and phenotype groups.
# Reads the consolidated p-value matrix produced by COMBINE_GWAS.
#
# Input:
#   gwas_pvals_<name>.csv  -- space-separated; columns:
#                             chr rs ps allele1 allele0
#                             p_lrt_<trait1> p_score_<trait1>
#                             p_lrt_<trait2> p_score_<trait2> ...
#
# Significance tiers (two optional thresholds):
#   sig_strict  -- p below the stricter threshold  → SIG_STRICT_COLOR
#   sig_lenient -- p below only the lenient threshold → SIG_LENIENT_COLOR
#   background  -- alternating grey by chromosome
#
# Outputs (one PDF per trait, one per group):
#   Manhattan_<name>_<papername>.pdf
#   Manhattan_<name>_group_<group>.pdf

# ── Sourcing ───────────────────────────────────────────────────────────────────
# LOCAL <- TRUE  : run interactively; argparse skipped, paths from args list below.
# LOCAL <- FALSE : run via Rscript on HPC / inside Nextflow.
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
    library(data.table)
    library(ggrepel)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
if (LOCAL) {
    args <- list(
        pvals          = "output/gwas/gwas_pvals_anxfb_test.csv",
        name           = "anxfb_test",
        yaml           = "data/anxfb.yaml",
        outdir         = "output/plots/manhattan",
        threshold_file = NULL,    # path to permutation threshold file, or NULL
        threshold      = 5e-8,    # fallback genome-wide threshold
        palette        = "nejm",
        annot_file     = "output/annotations/gene_annotations.csv",    # set to "output/annotations/gene_annotations.csv" to enable labels
        label_gap      = 1000     # minimum kb between labeled loci
    )
} else {
    parser <- ArgumentParser()
    parser$add_argument("--pvals", "-p", required = TRUE,
        help = "gwas_pvals_<name>.csv from COMBINE_GWAS")
    parser$add_argument("--name", "-n", required = TRUE,
        help = "Run name (params.name from Nextflow, e.g. anxfb_test)")
    parser$add_argument("--yaml", "-y", required = TRUE,
        help = "Pipeline YAML with phenotype groups and paper names")
    parser$add_argument("--outdir", "-o", default = ".",
        help = "Output directory (default: current dir)")
    parser$add_argument("--threshold_file",
        help = "Single-value file with permutation p-value threshold (from CALC_PERM_THRESHOLD)")
    parser$add_argument("--threshold", type = "double", default = 5e-8,
        help = "Fallback genome-wide significance threshold (default: 5e-8); used when --threshold_file is absent")
    parser$add_argument("--palette", default = "nejm",
        choices = c("nejm", "npg", "bmj", "jco", "lancet", "jama"),
        help = "ggsci palette for group colors in the combined plot legend (default: nejm)")
    parser$add_argument("--annot_file",
        help = "gene_annotations.csv from GET_GENES; enables MGI symbol labels on sig SNPs (single-trait plots only)")
    parser$add_argument("--label_gap", type = "double", default = 1000,
        help = "Min distance in kb between labeled loci (default: 1000)")
    args <- parser$parse_args()
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# ── Phenotype metadata ─────────────────────────────────────────────────────────
cfg              <- load_phenotype_config(args$yaml)
colors           <- assign_group_colors(cfg$pnames, cfg$groups_order, palette = args$palette)
pnames           <- colors$pnames
rownames(pnames) <- cfg$pheno_names

# ── Significance thresholds ────────────────────────────────────────────────────
# perm_thresh: from permutation testing (stricter, solid line)
# default_thresh: genome-wide fallback (dashed line)
# Either or both may be NULL/absent.
perm_thresh <- if (!is.null(args$threshold_file) &&
                   file.exists(args$threshold_file))
    suppressWarnings(as.numeric(readLines(args$threshold_file, warn = FALSE)[1])) else NULL

default_thresh <- args$threshold

# ── Gene annotation lookup (optional) ─────────────────────────────────────────
# Produces gene_labels: data.table(rs, mgi_symbol) — one best gene per SNP.
# SNP key: prefer lead_snp column (present in susie mode = one row per locus);
#          fall back to snp column (clump mode).
# Gene preference: protein_coding > other biotypes, then closest to SNP.
gene_labels <- NULL
if (!is.null(args$annot_file) && file.exists(args$annot_file)) {
    annot <- fread(args$annot_file)

    snp_col <- if ("lead_snp" %in% names(annot) && !all(is.na(annot$lead_snp)))
                   "lead_snp" else "snp"
    setnames(annot, snp_col, "rs_key")

    # Distance from SNP position to gene body (0 if SNP is inside the gene)
    annot[, dist_to_gene := ifelse(
        ps >= start_position & ps <= end_position, 0L,
        pmin(abs(ps - start_position), abs(ps - end_position))
    )]

    gene_labels <- annot[
        mgi_symbol != "" & !is.na(mgi_symbol)
    ][
        order(rs_key,
              -(gene_biotype == "protein_coding"),   # protein_coding first
              dist_to_gene)
    ][
        !duplicated(rs_key),
        .(rs = rs_key, mgi_symbol)
    ]
    message(sprintf("Gene labels loaded: %d annotated SNPs", nrow(gene_labels)))
}

# ── Load p-values ──────────────────────────────────────────────────────────────
pvals <- fread(args$pvals, sep = " ")
pvals[, chr := ifelse(chr == "X", 23L, suppressWarnings(as.integer(chr)))]
pvals <- pvals[!is.na(chr)]

p_lrt_cols  <- grep("^p_lrt_",   names(pvals), value = TRUE)
trait_names <- sub("^p_lrt_", "", p_lrt_cols)

# ── Cumulative genomic positions ───────────────────────────────────────────────
chr_offsets <- pvals[, .(chr_len = max(as.numeric(ps))), by = chr][order(chr)]
chr_offsets[, offset := cumsum(shift(chr_len, fill = 0))]
pvals <- merge(pvals, chr_offsets[, .(chr, offset)], by = "chr")
pvals[, ps_cum := as.numeric(ps) + offset]

axis_df <- pvals[, .(center = (min(ps_cum) + max(ps_cum)) / 2), by = chr][order(chr)]
axis_df[, label := ifelse(chr == 23L, "X", as.character(chr))]

# ── Core Manhattan plotting function ──────────────────────────────────────────
# dat        : data.table with columns rs (character), chr (integer),
#              ps_cum (numeric), p (numeric)
# title      : plot title string
# add_labels : if TRUE and gene_labels is loaded, annotate significant peaks
#              with MGI symbols via ggrepel (single-trait plots only)
make_manhattan <- function(dat, title, add_labels = FALSE) {
    dat <- dat[!is.na(p) & p > 0]

    # Thresholds sorted so thr[1] is the strictest (smallest p)
    thr <- sort(c(perm_thresh, default_thresh))

    dat[, sig_cat := ifelse(chr %% 2L == 0L, "chr_even", "chr_odd")]
    if (length(thr) >= 2L) dat[p <= thr[2L], sig_cat := "sig_lenient"]
    if (length(thr) >= 1L) dat[p <= thr[1L], sig_cat := "sig_strict"]
    dat[, sig_cat := factor(sig_cat,
                            levels = c("chr_odd", "chr_even",
                                       "sig_lenient", "sig_strict"))]

    color_map <- c(chr_odd     = CHROM_COLORS[1],
                   chr_even    = CHROM_COLORS[2],
                   sig_lenient = SIG_LENIENT_COLOR,
                   sig_strict  = SIG_STRICT_COLOR)

    size_map  <- c(chr_odd     = 0.6,
                   chr_even    = 0.6,
                   sig_lenient = 1.5,
                   sig_strict  = 1.5)

    plt <- ggplot(dat, aes(x = ps_cum, y = -log10(p),
                           color = sig_cat, size = sig_cat)) +
        geom_point(alpha = 0.85) +
        scale_color_manual(values = color_map, drop = FALSE) +
        scale_size_manual(values  = size_map,  drop = FALSE) +
        scale_x_continuous(breaks = axis_df$center,
                           labels = axis_df$label,
                           expand = c(0.01, 0)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(x     = "Chromosome",
             y     = expression(-log[10](italic(p))),
             title = title) +
        theme_manhattan()

    if (!is.null(perm_thresh))
        plt <- plt + geom_hline(yintercept = -log10(perm_thresh),
                                color = "black", linetype = "solid",
                                linewidth = 0.5)
    if (!is.null(default_thresh))
        plt <- plt + geom_hline(yintercept = -log10(default_thresh),
                                color = "black", linetype = "dashed",
                                linewidth = 0.5)

    # ── Gene labels on significant peaks ──────────────────────────────────────
    if (add_labels && !is.null(gene_labels) && length(thr) >= 1L) {

        # Significant SNPs that have a gene annotation
        sig_dat <- merge(dat[p <= thr[1L] & sig_cat == "sig_strict"],
                         gene_labels, by = "rs", all.x = FALSE)

        if (nrow(sig_dat) > 0L) {
            # Greedy peak selection: sort by p, skip any SNP within label_gap kb
            # of an already-selected peak (uses cumulative position, safe across chrs)
            sig_dat  <- sig_dat[order(p)]
            keep     <- logical(nrow(sig_dat))
            chosen_x <- numeric(0)
            gap_bp   <- args$label_gap * 1000

            for (i in seq_len(nrow(sig_dat))) {
                if (!length(chosen_x) ||
                    all(abs(sig_dat$ps_cum[i] - chosen_x) > gap_bp)) {
                    keep[i]  <- TRUE
                    chosen_x <- c(chosen_x, sig_dat$ps_cum[i])
                }
            }
            label_dat <- sig_dat[keep]

            plt <- plt +
                ggrepel::geom_label_repel(
                    data             = label_dat,
                    aes(x = ps_cum, y = -log10(p), label = mgi_symbol),
                    inherit.aes      = FALSE,
                    size             = 2.5,
                    box.padding      = 0.4,
                    point.padding    = 0.3,
                    label.size       = 0.2,
                    fill             = "white",
                    alpha            = 0.85,
                    colour           = "grey20",
                    segment.color    = "grey50",
                    min.segment.length = 0.1,
                    max.overlaps     = Inf
                )
        }
    }

    plt
}

# ── Winner-colored group Manhattan (significant SNPs colored by winning trait) ──
# dat_full   : data.table with rs, chr, ps_cum + p_lrt_<trait> columns for group
# grp_traits : character vector of pheno_names keys (row names in pnames)
# title      : plot title string
make_manhattan_winner <- function(dat_full, grp_traits, title) {
    trait_cols        <- intersect(paste0("p_lrt_", grp_traits), names(dat_full))
    if (length(trait_cols) == 0L) return(invisible(NULL))
    grp_traits_active <- sub("^p_lrt_", "", trait_cols)

    paper_names <- pnames[grp_traits_active, "PaperName"]
    names(paper_names) <- grp_traits_active

    # Per-SNP minimum p and winning PaperName
    p_mat       <- as.matrix(dat_full[, trait_cols, with = FALSE])
    min_p       <- apply(p_mat, 1L, min, na.rm = TRUE)
    win_col_idx <- apply(p_mat, 1L, function(x)
                         if (all(is.na(x))) 1L else which.min(x))
    win_paper   <- paper_names[grp_traits_active[win_col_idx]]

    thr <- sort(c(perm_thresh, default_thresh))

    dat <- dat_full[, .(rs, chr, ps_cum)]
    dat[, min_p     := min_p]
    dat[, win_paper := win_paper]
    dat[, sig_cat   := ifelse(chr %% 2L == 0L, "chr_even", "chr_odd")]
    if (length(thr) >= 2L) dat[min_p <= thr[2L], sig_cat := "sig_lenient"]
    if (length(thr) >= 1L) dat[min_p <= thr[1L], sig_cat := win_paper]

    unique_papers <- unique(paper_names)
    active_papers <- intersect(unique_papers, unique(dat$sig_cat))

    n_traits     <- length(unique_papers)
    trait_pal    <- cluster_colors(n_traits, palette = if (n_traits <= 10L) "d3" else "igv")
    trait_colors <- setNames(trait_pal, unique_papers)

    all_levels <- c("chr_odd", "chr_even", "sig_lenient", unique_papers)
    dat[, sig_cat := factor(sig_cat, levels = all_levels)]

    color_map <- c(chr_odd     = CHROM_COLORS[1],
                   chr_even    = CHROM_COLORS[2],
                   sig_lenient = SIG_LENIENT_COLOR,
                   trait_colors)

    size_map <- c(chr_odd     = 0.6,
                  chr_even    = 0.6,
                  sig_lenient = 1.5,
                  setNames(rep(1.5, n_traits), unique_papers))

    plt <- ggplot(dat[!is.na(min_p) & min_p > 0],
                  aes(x = ps_cum, y = -log10(min_p),
                      color = sig_cat, size = sig_cat)) +
        geom_point(alpha = 0.85) +
        scale_color_manual(
            values = color_map,
            breaks = active_papers,
            drop   = FALSE
        ) +
        scale_size_manual(
            values = size_map,
            breaks = active_papers,
            drop   = FALSE
        ) +
        scale_x_continuous(breaks = axis_df$center,
                           labels = axis_df$label,
                           expand = c(0.01, 0)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(x     = "Chromosome",
             y     = expression(-log[10](italic(p))),
             color = "Trait",
             size  = "Trait",
             title = title) +
        theme_manhattan() +
        theme(legend.position = if (length(active_papers) > 0L) "right" else "none")

    if (!is.null(perm_thresh))
        plt <- plt + geom_hline(yintercept = -log10(perm_thresh),
                                color = "black", linetype = "solid",
                                linewidth = 0.5)
    if (!is.null(default_thresh))
        plt <- plt + geom_hline(yintercept = -log10(default_thresh),
                                color = "black", linetype = "dashed",
                                linewidth = 0.5)
    plt
}

# ── Combined all-groups Manhattan (significant SNPs colored by winning group) ──
make_manhattan_combined <- function() {
    if (length(p_lrt_cols) == 0L) return(invisible(NULL))

    grpcol_map <- colors$grpcol   # named vector: group -> hex color

    # Per-SNP minimum p across all traits and the winning group
    p_mat       <- as.matrix(pvals[, p_lrt_cols, with = FALSE])
    min_p       <- apply(p_mat, 1L, min, na.rm = TRUE)
    win_col_idx <- apply(p_mat, 1L, function(x)
                         if (all(is.na(x))) 1L else which.min(x))
    win_trait   <- trait_names[win_col_idx]
    win_group   <- pnames[win_trait, "Group"]

    thr <- sort(c(perm_thresh, default_thresh))

    dat <- pvals[, .(rs, chr, ps_cum)]
    dat[, min_p     := min_p]
    dat[, win_group := win_group]
    dat[, sig_cat   := ifelse(chr %% 2L == 0L, "chr_even", "chr_odd")]
    if (length(thr) >= 2L) dat[min_p <= thr[2L], sig_cat := "sig_lenient"]
    if (length(thr) >= 1L) dat[min_p <= thr[1L], sig_cat := win_group]

    active_groups <- intersect(names(grpcol_map), unique(dat$sig_cat))

    all_levels <- c("chr_odd", "chr_even", "sig_lenient", names(grpcol_map))
    dat[, sig_cat := factor(sig_cat, levels = all_levels)]

    color_map <- c(chr_odd     = CHROM_COLORS[1],
                   chr_even    = CHROM_COLORS[2],
                   sig_lenient = SIG_LENIENT_COLOR,
                   grpcol_map)

    size_map <- c(chr_odd     = 0.6,
                  chr_even    = 0.6,
                  sig_lenient = 1.5,
                  setNames(rep(1.5, length(grpcol_map)), names(grpcol_map)))

    plt <- ggplot(dat[!is.na(min_p) & min_p > 0],
                  aes(x = ps_cum, y = -log10(min_p),
                      color = sig_cat, size = sig_cat)) +
        geom_point(alpha = 0.85) +
        scale_color_manual(
            values = color_map,
            breaks = active_groups,
            drop   = FALSE
        ) +
        scale_size_manual(
            values = size_map,
            breaks = active_groups,
            drop   = FALSE
        ) +
        scale_x_continuous(breaks = axis_df$center,
                           labels = axis_df$label,
                           expand = c(0.01, 0)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(x     = "Chromosome",
             y     = expression(-log[10](italic(p))),
             color = "Group",
             size  = "Group",
             title = paste0(args$name, " — all traits")) +
        theme_manhattan() +
        theme(legend.position = if (length(active_groups) > 0L) "right" else "none")

    if (!is.null(perm_thresh))
        plt <- plt + geom_hline(yintercept = -log10(perm_thresh),
                                color = "black", linetype = "solid",
                                linewidth = 0.5)
    if (!is.null(default_thresh))
        plt <- plt + geom_hline(yintercept = -log10(default_thresh),
                                color = "black", linetype = "dashed",
                                linewidth = 0.5)
    plt
}

# ── Individual trait Manhattan plots ──────────────────────────────────────────
message("Plotting ", length(p_lrt_cols), " individual Manhattan plots...")

for (i in seq_along(p_lrt_cols)) {
    trait     <- trait_names[i]
    papername <- pnames[trait, "PaperName"]

    if (is.na(papername)) {
        message("  Skipping ", trait, " (not found in YAML)")
        next
    }

    dat <- pvals[, .(rs, chr, ps_cum, p = get(p_lrt_cols[i]))]
    plt <- make_manhattan(dat, title = papername, add_labels = !is.null(gene_labels))

    save_plot(plt,
              path = file.path(args$outdir,
                               paste0("Manhattan_", args$name, "_", papername, ".pdf")))
    message("  ", papername)
}

# ── Group Manhattan plots (winner-colored by trait within group) ───────────────
all_groups <- unique(pnames$Group)
message("Plotting ", length(all_groups), " group Manhattan plots (winner-colored by trait)...")

for (grp in all_groups) {
    grp_traits <- rownames(pnames)[pnames$Group == grp]
    grp_cols   <- intersect(paste0("p_lrt_", grp_traits), names(pvals))

    if (length(grp_cols) == 0) {
        message("  Skipping group '", grp, "' (no matching p_lrt columns)")
        next
    }

    dat_full <- pvals[, c("rs", "chr", "ps_cum", grp_cols), with = FALSE]
    plt      <- make_manhattan_winner(dat_full, grp_traits, title = grp)
    out_stem <- gsub("[^A-Za-z0-9]+", "_", grp)
    save_plot(plt,
              path  = file.path(args$outdir,
                                paste0("Manhattan_", args$name, "_group_", out_stem, ".pdf")),
              width = FIG$full_width + 1.5)
    message("  ", grp)
}

# ── Combined all-groups Manhattan (winner-colored by group) ───────────────────
message("Plotting combined all-groups Manhattan (winner-colored by group)...")
plt_combined <- make_manhattan_combined()
if (!is.null(plt_combined)) {
    save_plot(plt_combined,
              path  = file.path(args$outdir,
                                paste0("Manhattan_", args$name, "_combined.pdf")),
              width = FIG$full_width + 1.5)
    message("  combined")
}

message("Done. Plots written to: ", normalizePath(args$outdir))
