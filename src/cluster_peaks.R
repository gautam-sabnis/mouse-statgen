#!/usr/bin/env Rscript
#
# cluster_peaks.R
#
# Clusters fine-mapped GWAS loci by cross-trait p-value patterns (Components 6 + 7).
#
# Loci are defined as the lead SNPs from PLINK clumping that produced at least
# one SuSiE credible set. A pleiotropy matrix (loci x traits, -log10 p-values)
# is built from gwas_pvals, clustered by k-means, and visualised as a heatmap.
# Cluster-colored Manhattan plots are then produced for each trait and group.
#
# Inputs:
#   gwas_pvals_<name>.csv          from COMBINE_GWAS
#   trait*_locus_*.txt files       from SUSIE_FINEMAP (passed as a list)
#   Pipeline YAML
#
# Outputs:
#   pleiotropy_heatmap_<name>.pdf/.svg     loci x traits heatmap, rows by cluster
#   cluster_assignments_<name>.csv         lead_snp, chr, bp, cluster
#   Manhattan_clusters_<name>_<trait>.pdf  per-trait, lead SNPs colored by cluster
#   Manhattan_clusters_<name>_group_<grp>.pdf  per-group, lead SNPs colored by cluster

# ── Sourcing ───────────────────────────────────────────────────────────────────
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
    library(gplots)
    library(grid)
    library(UpSetR)
    library(ggVennDiagram)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
if (LOCAL) {
    args <- list(
        pvals       = "output/gwas/gwas_pvals_anxfb_test.csv",
        susie       = list.files("output/susie", pattern = "trait.*_locus_.*.txt",
                                 full.names = TRUE),
        name        = "anxfb_test",
        yaml        = "data/anxfb.yaml",
        outdir      = "output/plots/clusters",
        clusters    = 7L,
        pvalthr     = 5,        # -log10(p) threshold for significance line on Manhattans
        palette     = "nejm",
        clumped_dir = "output/clump_combined",  # Part A: per-group CLUMP counts + overlap
        clump_dir   = "output/clump",           # Part D: per-trait CLUMP counts
        overlap_kb  = 1000                      # Part E2: proximity window for CLUMP overlap
    )
} else {
    parser <- ArgumentParser()
    parser$add_argument("--pvals", "-p", required = TRUE,
        help = "gwas_pvals_<name>.csv from COMBINE_GWAS")
    parser$add_argument("--susie", "-s", required = TRUE, nargs = "+",
        help = "SuSiE locus files (trait*_locus_*.txt) from SUSIE_FINEMAP")
    parser$add_argument("--name", "-n", required = TRUE,
        help = "Run name (params.name from Nextflow, e.g. anxfb_test)")
    parser$add_argument("--yaml", "-y", required = TRUE,
        help = "Pipeline YAML with phenotype groups and paper names")
    parser$add_argument("--outdir", "-o", default = ".",
        help = "Output directory (default: current dir)")
    parser$add_argument("--clusters", "-k", type = "integer", default = 7L,
        help = "Number of k-means clusters (default: 7)")
    parser$add_argument("--pvalthr", type = "double", default = 5,
        help = "-log10(p) threshold for significance line on cluster Manhattans (default: 5)")
    parser$add_argument("--palette", default = "nejm",
        choices = c("nejm", "npg", "bmj", "jco", "lancet", "jama"),
        help = "ggsci palette for group column side colors (default: nejm)")
    parser$add_argument("--clumped_dir", default = NULL,
        help = "Directory with .clumped files from CLUMP_COMBINED (enables per-group QTL count bar chart)")
    parser$add_argument("--clump_dir", default = NULL,
        help = "Directory with trait<N>.clumped files from CLUMP (enables per-trait CLUMP QTL count bar chart)")
    parser$add_argument("--overlap_kb", type = "double", default = 1000,
        help = "Proximity window in kb for CLUMP-based QTL overlap between groups (default: 1000)")
    args <- parser$parse_args()
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(490)   # match original for reproducibility

# ── Phenotype metadata ─────────────────────────────────────────────────────────
cfg              <- load_phenotype_config(args$yaml)
colors           <- assign_group_colors(cfg$pnames, cfg$groups_order,
                                        palette = args$palette)
pnames           <- colors$pnames
rownames(pnames) <- cfg$pheno_names
grpcol           <- colors$grpcol

# ── Load p-value matrix ────────────────────────────────────────────────────────
pvals      <- fread(args$pvals, sep = " ")
p_lrt_cols <- grep("^p_lrt_", names(pvals), value = TRUE)
trait_names <- sub("^p_lrt_", "", p_lrt_cols)

# ── Collect SuSiE credible set lead SNPs ──────────────────────────────────────
# A lead SNP qualifies if it produced at least one SNP with cs != NA.
susie_all <- rbindlist(lapply(args$susie, fread), fill = TRUE)

# Lead SNPs with at least one credible set member
cs_leads <- susie_all[!is.na(cs), .(chr = chr[1], bp = bp[1]), by = lead_snp]

if (nrow(cs_leads) == 0)
    stop("No SuSiE credible sets found in the provided locus files.")

message(sprintf("%d loci with credible sets across all traits", nrow(cs_leads)))

# ── Build pleiotropy matrix (loci x traits, -log10 p_lrt) ─────────────────────
# Filter p-value matrix to lead SNPs present in gwas_pvals
pvals_leads <- pvals[rs %in% cs_leads$lead_snp,
                     c("rs", p_lrt_cols), with = FALSE]

# Some lead SNPs may be absent from gwas_pvals (e.g. multi-allelic or filtered)
missing <- setdiff(cs_leads$lead_snp, pvals_leads$rs)
if (length(missing) > 0)
    message(sprintf("  %d lead SNP(s) not found in gwas_pvals and will be dropped: %s",
                    length(missing), paste(missing, collapse = ", ")))

pmat <- as.matrix(pvals_leads[, p_lrt_cols, with = FALSE])
rownames(pmat) <- pvals_leads$rs

# Transform to -log10; replace NA/Inf with 0 (no evidence)
pmat <- -log10(pmat)
pmat[!is.finite(pmat)] <- 0

# ── k-means clustering ─────────────────────────────────────────────────────────
kk <- kmeans(pmat, centers = args$clusters, nstart = 25)

# Within each cluster, order rows by hierarchical clustering
rowarr <- character(0)
for (cl in seq_len(args$clusters)) {
    idx <- which(kk$cluster == cl)
    if (length(idx) == 1) {
        rowarr <- c(rowarr, rownames(pmat)[idx])
    } else {
        hc     <- hclust(dist(pmat[idx, , drop = FALSE]))
        rowarr <- c(rowarr, rownames(pmat)[idx][hc$order])
    }
}
pmat_ord <- pmat[rowarr, , drop = FALSE]

# ── Cluster and group color bars ───────────────────────────────────────────────
ccols    <- cluster_colors(args$clusters)
row_cols <- ccols[kk$cluster[rowarr]]   # one color per row (cluster)

# Column side colors: group color for each trait (in p_lrt column order)
col_cols <- vapply(trait_names, function(tr) {
    col <- pnames[tr, "color"]
    if (is.na(col)) "#CCCCCC" else col
}, character(1))

# ── Heatmap ────────────────────────────────────────────────────────────────────
hwid <- if (ncol(pmat_ord) > 40) FIG$full_width * 2 else FIG$full_width

plot_heatmap <- function() {
    heatmap.2(
        pmat_ord,
        col          = HEATMAP_COLORS,
        Rowv         = FALSE,          # rows already ordered
        Colv         = TRUE,           # cluster columns
        dendrogram   = "col",
        scale        = "none",
        trace        = "none",
        RowSideColors = row_cols,
        ColSideColors = col_cols,
        labRow       = NA,             # too many rows to label
        hclustfun    = function(x) hclust(x, method = "average"),
        distfun      = function(x) dist(scale(x)),
        margins      = c(12, 8),
        srtCol       = 45,
        key          = TRUE,
        key.title    = "-log10(p)",
        density.info = "none"
    )
}

out_pdf <- file.path(args$outdir, paste0("pleiotropy_heatmap_", args$name, ".pdf"))
out_svg <- file.path(args$outdir, paste0("pleiotropy_heatmap_", args$name, ".svg"))

PLOT_DEVICE(out_pdf, width = hwid, height = FIG$panel_height + 1, family = FIG_FONT)
hplt <- plot_heatmap()
dev.off()

svg(out_svg, width = hwid, height = FIG$panel_height + 1, family = FIG_FONT)
eval(hplt$call)
dev.off()

message("Written: ", out_pdf)

# ── Save cluster assignments ───────────────────────────────────────────────────
cluster_tbl <- cs_leads[lead_snp %in% rownames(pmat_ord)]
cluster_tbl[, cluster := kk$cluster[rowarr][match(lead_snp, rowarr)]]

out_csv <- file.path(args$outdir, paste0("cluster_assignments_", args$name, ".csv"))
fwrite(cluster_tbl, out_csv)
message("Written: ", out_csv)

# ── Component 7: Cluster-colored Manhattan plots ───────────────────────────────
# Cumulative genomic positions (same logic as plot_manhattan.R)
pvals[, chr_int := ifelse(chr == "X", 23L, suppressWarnings(as.integer(chr)))]
pvals <- pvals[!is.na(chr_int)]

chr_offsets <- pvals[, .(chr_len = max(as.numeric(ps))), by = chr_int][order(chr_int)]
chr_offsets[, offset := cumsum(shift(chr_len, fill = 0))]
pvals <- merge(pvals, chr_offsets[, .(chr_int, offset)], by = "chr_int")
pvals[, ps_cum := as.numeric(ps) + offset]

axis_df <- pvals[, .(center = (min(ps_cum) + max(ps_cum)) / 2), by = chr_int][order(chr_int)]
axis_df[, label := ifelse(chr_int == 23L, "X", as.character(chr_int))]

# Attach cluster assignment to every lead SNP row in pvals
pvals <- merge(pvals,
               cluster_tbl[, .(rs = lead_snp, cluster = as.factor(cluster))],
               by = "rs", all.x = TRUE)

make_cluster_manhattan <- function(dat, title) {
    dat <- dat[!is.na(p) & p > 0]

    # Points: lead SNPs get cluster color; all others alternating grey
    dat[, point_col := ifelse(chr_int %% 2L == 0L, CHROM_COLORS[2], CHROM_COLORS[1])]
    dat[!is.na(cluster), point_col := ccols[as.integer(cluster)]]
    dat[, point_size := ifelse(!is.na(cluster), 1.8, 0.6)]

    ggplot(dat, aes(x = ps_cum, y = -log10(p))) +
        geom_point(color = dat$point_col, size = dat$point_size, alpha = 0.85) +
        geom_hline(yintercept = args$pvalthr,
                   color = SIG_LINE_COLOR, linetype = "solid", linewidth = 0.5) +
        scale_x_continuous(breaks = axis_df$center,
                           labels = axis_df$label,
                           expand = c(0.01, 0)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(x = "Chromosome",
             y = expression(-log[10](italic(p))),
             title = title) +
        theme_manhattan()
}

# Per-trait cluster-colored Manhattans
message("Plotting ", length(p_lrt_cols), " cluster-colored individual Manhattan plots...")

for (i in seq_along(p_lrt_cols)) {
    trait     <- trait_names[i]
    papername <- pnames[trait, "PaperName"]
    if (is.na(papername)) next

    dat <- pvals[, .(rs, chr_int, ps_cum, p = get(p_lrt_cols[i]), cluster)]
    plt <- make_cluster_manhattan(dat, title = papername)
    save_plot(plt,
              path = file.path(args$outdir,
                               paste0("Manhattan_clusters_", args$name,
                                      "_", papername, ".pdf")))
    message("  ", papername)
}

# Per-group cluster-colored Manhattans (min p-value across group traits)
all_groups <- unique(pnames$Group)
message("Plotting ", length(all_groups), " cluster-colored group Manhattan plots...")

for (grp in all_groups) {
    grp_traits <- rownames(pnames)[pnames$Group == grp]
    grp_cols   <- intersect(paste0("p_lrt_", grp_traits), names(pvals))
    if (length(grp_cols) == 0) next

    grp_mat  <- as.matrix(pvals[, grp_cols, with = FALSE])
    min_pval <- apply(grp_mat, 1L, min, na.rm = TRUE)
    dat      <- pvals[, .(rs, chr_int, ps_cum, p = min_pval, cluster)]

    plt      <- make_cluster_manhattan(dat, title = grp)
    out_stem <- gsub("[^A-Za-z0-9]+", "_", grp)
    save_plot(plt,
              path = file.path(args$outdir,
                               paste0("Manhattan_clusters_", args$name,
                                      "_group_", out_stem, ".pdf")))
    message("  ", grp)
}

message("Done. All outputs written to: ", normalizePath(args$outdir))

# ── Component 8: Pleiotropy analysis + QTL counts ─────────────────────────────
#
# Part B (always): For each SuSiE lead SNP, count which phenotype groups have
#   at least one trait with p < threshold. → pleiotropy_summary CSV + histogram.
# Part C (always): Count SuSiE credible-set loci per trait and per group.
#   → qtl_counts_susie_<name>.pdf
# Part A (requires --clumped_dir): Per-group QTL counts from CLUMP_COMBINED.
#   → qtl_counts_groups_<name>.pdf
# Part D (requires --clump_dir): Per-trait QTL counts from CLUMP lead SNPs.
#   → qtl_counts_clump_<name>.pdf

message("── Pleiotropy analysis + QTL counts ────────────────────────────────")

pthresh_p <- 10^(-args$pvalthr)   # convert -log10 threshold back to p-value

# ── Shared helpers ─────────────────────────────────────────────────────────────

# Read a PLINK .clumped file; returns 0-row data.table if absent/empty/no SNPs
read_clumped <- function(path) {
    if (!file.exists(path) || file.info(path)$size == 0L) return(data.table())
    dt <- fread(path, fill = TRUE)
    if (nrow(dt) == 0L || !"SNP" %in% names(dt)) return(data.table())
    dt[, .(SNP, CHR, BP, P)]
}

# Trait metadata: one row per trait in YAML order.
# group_rank preserves cfg$groups_order for consistent sort across all charts.
group_rank_vec <- setNames(seq_along(cfg$groups_order), cfg$groups_order)
trait_map <- data.table(
    trait     = seq_along(cfg$pheno_names),
    PaperName = pnames[cfg$pheno_names, "PaperName"],
    Group     = pnames[cfg$pheno_names, "Group"]
)[, grp_rank := group_rank_vec[Group]]

# Plot height scaled to number of bars (~0.22 in per bar, minimum one panel)
pvh_trait <- max(FIG$panel_height, 0.22 * nrow(trait_map))

# Horizontal bar chart: group-colored, coord_flip, count labels
bar_qtl <- function(dt, title) {
    # dt must have: PaperName (ordered factor), n_qtl (integer), Group (character)
    ggplot(dt, aes(x = PaperName, y = n_qtl, fill = Group)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
        geom_text(aes(label = n_qtl), hjust = -0.15, size = 2.3) +
        scale_fill_manual(values = grpcol) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
        coord_flip() +
        labs(x = NULL, y = "QTLs", title = title) +
        theme_pve() +
        theme(legend.position = "none")
}

# ── Part B: Group-level significance at SuSiE lead SNPs ──────────────────────
# pvals already has chr_int, ps_cum, and cluster columns from the Manhattan section.
leads_pv <- pvals[rs %in% rownames(pmat_ord)]

sig_mat <- matrix(FALSE,
                  nrow     = nrow(leads_pv),
                  ncol     = length(all_groups),
                  dimnames = list(leads_pv$rs, all_groups))

for (grp in all_groups) {
    grp_traits <- rownames(pnames)[pnames$Group == grp]
    grp_cols   <- intersect(paste0("p_lrt_", grp_traits), names(leads_pv))
    if (length(grp_cols) == 0) next
    grp_mat  <- as.matrix(leads_pv[, grp_cols, with = FALSE])
    min_pval <- apply(grp_mat, 1L, min, na.rm = TRUE)
    sig_mat[leads_pv$rs, grp] <- min_pval < pthresh_p
}

plei_tbl <- data.table(
    lead_snp     = rownames(sig_mat),
    n_groups     = rowSums(sig_mat, na.rm = TRUE),
    which_groups = apply(sig_mat, 1L, function(x) paste(names(which(x)), collapse = "|"))
)
plei_tbl <- merge(plei_tbl, cs_leads[, .(lead_snp, chr, bp)],
                  by = "lead_snp", all.x = TRUE)
plei_tbl <- merge(plei_tbl, cluster_tbl[, .(lead_snp, cluster)],
                  by = "lead_snp", all.x = TRUE)
plei_tbl <- plei_tbl[order(-n_groups, chr, bp)]

out_plei <- file.path(args$outdir, paste0("pleiotropy_summary_", args$name, ".csv"))
fwrite(plei_tbl, out_plei)
message("Written: ", out_plei)

# Histogram: QTLs by number of groups they are significant for
hist_plt <- ggplot(plei_tbl, aes(x = n_groups, fill = as.factor(cluster))) +
    geom_bar(width = 0.7, color = "white") +
    scale_fill_manual(values = ccols, name = "Cluster",
                      breaks = seq_len(args$clusters)) +
    scale_x_continuous(breaks = 0:length(all_groups), minor_breaks = NULL) +
    labs(x     = sprintf("Groups with p < 10^%g", -args$pvalthr),
         y     = "QTLs",
         title = "Pleiotropy: group count per QTL") +
    theme_gwas()

save_plot(hist_plt,
          path   = file.path(args$outdir,
                             paste0("pleiotropy_hist_", args$name, ".pdf")),
          width  = FIG$half_width,
          height = FIG$panel_height)
message("Written: pleiotropy_hist_", args$name, ".pdf")

# ── Part C: Per-trait SuSiE credible-set locus counts (always) ───────────────
# Count unique lead SNPs with at least one credible set per trait.
# Traits with zero credible sets (no loci fine-mapped) are included as 0.
trait_susie <- merge(
    trait_map,
    susie_all[!is.na(cs), .(n_qtl = uniqueN(lead_snp)), by = trait],
    by = "trait", all.x = TRUE
)
trait_susie[is.na(n_qtl), n_qtl := 0L]

# Order: group position in groups_order, then n_qtl descending within group
trait_susie <- trait_susie[order(grp_rank, -n_qtl)]
trait_susie[, PaperName := factor(PaperName, levels = rev(PaperName))]

out_susie_trait <- file.path(args$outdir,
                             paste0("qtl_counts_susie_", args$name, ".pdf"))
save_plot(bar_qtl(trait_susie, "QTL counts per trait (SuSiE credible sets)"),
          path = out_susie_trait, width = FIG$half_width * 1.5, height = pvh_trait)
message("Written: ", out_susie_trait)

# Also save per-group SuSiE counts (sum over traits within each group)
grp_susie <- trait_susie[, .(n_qtl = sum(n_qtl), Group = Group[1L]), by = .(grp_rank)]
grp_susie[, PaperName := cfg$groups_order[grp_rank]]
grp_susie[, PaperName := factor(PaperName, levels = rev(cfg$groups_order))]

pvh_grp_s <- max(FIG$panel_height, 0.22 * nrow(grp_susie))
out_susie_grp <- file.path(args$outdir,
                           paste0("qtl_counts_susie_groups_", args$name, ".pdf"))
save_plot(bar_qtl(grp_susie, "QTL counts by group (SuSiE credible sets)"),
          path = out_susie_grp, width = FIG$half_width * 1.5, height = pvh_grp_s)
message("Written: ", out_susie_grp)

# ── Part A: Per-group QTL counts from CLUMP_COMBINED (optional) ──────────────
if (!is.null(args$clumped_dir)) {
    n_all <- nrow(read_clumped(file.path(args$clumped_dir, "all.clumped")))
    message(sprintf("All-traits combined CLUMP QTLs: %d", n_all))

    grp_clump <- rbindlist(lapply(cfg$groups_order, function(grp) {
        stem  <- gsub("[^A-Za-z0-9]+", "_", grp)
        fpath <- file.path(args$clumped_dir, paste0(stem, ".clumped"))
        data.table(PaperName = grp, Group = grp,
                   n_qtl = nrow(read_clumped(fpath)))
    }))
    grp_clump[, PaperName := factor(PaperName, levels = rev(cfg$groups_order))]

    pvh_grp_c <- max(FIG$panel_height, 0.22 * nrow(grp_clump))
    out_grp_clump <- file.path(args$outdir,
                               paste0("qtl_counts_groups_", args$name, ".pdf"))
    save_plot(
        bar_qtl(grp_clump,
                sprintf("QTL counts by group  (all-traits: %d)", n_all)),
        path = out_grp_clump, width = FIG$half_width * 1.5, height = pvh_grp_c
    )
    message("Written: ", out_grp_clump)
}

# ── Part D: Per-trait QTL counts from CLUMP (optional) ───────────────────────
if (!is.null(args$clump_dir)) {
    n_qtl_v <- vapply(trait_map$trait, function(t) {
        nrow(read_clumped(file.path(args$clump_dir,
                                    sprintf("trait%d.clumped", t))))
    }, integer(1))

    trait_clump <- copy(trait_map)[, n_qtl := n_qtl_v]
    trait_clump <- trait_clump[order(grp_rank, -n_qtl)]
    trait_clump[, PaperName := factor(PaperName, levels = rev(PaperName))]

    out_clump_trait <- file.path(args$outdir,
                                 paste0("qtl_counts_clump_", args$name, ".pdf"))
    save_plot(bar_qtl(trait_clump, "QTL counts per trait (PLINK clump lead SNPs)"),
              path = out_clump_trait, width = FIG$half_width * 1.5, height = pvh_trait)
    message("Written: ", out_clump_trait)
}

# ── Part E: QTL overlap between groups ────────────────────────────────────────
#
# E1 (always): UpSet + Venn from the SuSiE-based pleiotropy_summary already
#   computed in Part B.  Each set = group; elements = SuSiE lead SNPs where
#   the group's min p_lrt < threshold.
#   → qtl_overlap_upset_susie_<name>.pdf
#   → qtl_overlap_venn_susie_<name>.pdf   (only when n_active_groups ≤ 4)
#
# E2 (requires --clumped_dir): UpSet + Venn based on proximity of independent
#   CLUMP lead SNPs across groups.  A QTL from all.clumped is assigned to a
#   group if that group has a lead SNP within --overlap_kb on the same chr.
#   → qtl_overlap_upset_clump_<name>.pdf
#   → qtl_overlap_venn_clump_<name>.pdf   (only when n_active_groups ≤ 4)

# Shared saving wrapper for UpSetR (grid graphics, not ggplot).
# suppressWarnings() silences UpSetR's internal use of deprecated ggplot2
# aesthetics (aes_string, size for lines/element_line) — not our code.
save_upset <- function(upset_obj, path, width = FIG$full_width,
                       height = FIG$panel_height * 1.5) {
    PLOT_DEVICE(path, width = width, height = height, family = FIG_FONT)
    suppressWarnings(print(upset_obj))
    dev.off()
}

# ── Part E1: SuSiE-based overlap ──────────────────────────────────────────────
# Build one set per group from plei_tbl (computed in Part B).
susie_sets <- Filter(length, setNames(
    lapply(all_groups, function(g)
        plei_tbl[n_groups > 0 & grepl(g, which_groups, fixed = TRUE), lead_snp]),
    all_groups
))

if (length(susie_sets) >= 2) {
    susie_ud <- UpSetR::fromList(susie_sets)

    save_upset(
        suppressWarnings(UpSetR::upset(
            susie_ud,
            sets            = rev(names(susie_sets)),
            order.by        = "freq",
            mb.ratio        = c(0.65, 0.35),
            text.scale      = c(1.2, 1.1, 1, 1, 1.1, 0.9),
            point.size      = 2.8,
            line.size       = 0.8,
            mainbar.y.label = "Shared QTLs",
            sets.x.label    = "QTLs per group"
        )),
        path = file.path(args$outdir,
                         paste0("qtl_overlap_upset_susie_", args$name, ".pdf"))
    )
    message("Written: qtl_overlap_upset_susie_", args$name, ".pdf")

    if (length(susie_sets) <= 4) {
        venn_susie <- ggVennDiagram::ggVennDiagram(susie_sets, label_alpha = 0) +
            scale_fill_gradient(low = "white", high = SIG_STRICT_COLOR) +
            scale_color_manual(values = rep("grey40", length(susie_sets))) +
            theme(legend.position = "none")
        save_plot(venn_susie,
                  path   = file.path(args$outdir,
                                     paste0("qtl_overlap_venn_susie_", args$name, ".pdf")),
                  width  = FIG$half_width * 1.5,
                  height = FIG$panel_height * 1.2)
        message("Written: qtl_overlap_venn_susie_", args$name, ".pdf")
    }
} else {
    message("Fewer than 2 groups with SuSiE QTLs — skipping E1 overlap plots")
}

# ── Part E2: CLUMP proximity-based overlap (optional) ─────────────────────────
if (!is.null(args$clumped_dir)) {
    all_qtl <- read_clumped(file.path(args$clumped_dir, "all.clumped"))

    if (nrow(all_qtl) >= 2) {
        overlap_bp <- args$overlap_kb * 1000

        # Per-group lead SNP positions (already read for Part A, re-read here)
        grp_pos <- setNames(
            lapply(all_groups, function(g) {
                stem <- gsub("[^A-Za-z0-9]+", "_", g)
                read_clumped(file.path(args$clumped_dir, paste0(stem, ".clumped")))
            }),
            all_groups
        )

        # For each all-traits QTL, flag which groups have a nearby lead SNP
        mem_mat <- matrix(FALSE,
                          nrow     = nrow(all_qtl),
                          ncol     = length(all_groups),
                          dimnames = list(all_qtl$SNP, all_groups))

        for (i in seq_len(nrow(all_qtl))) {
            for (g in all_groups) {
                gdt <- grp_pos[[g]]
                if (nrow(gdt) == 0) next
                mem_mat[i, g] <- any(gdt$CHR == all_qtl$CHR[i] &
                                     abs(gdt$BP - all_qtl$BP[i]) <= overlap_bp)
            }
        }

        # Named list of SNP IDs per group (for UpSetR / ggVennDiagram)
        clump_sets <- Filter(length, setNames(
            lapply(all_groups, function(g) all_qtl$SNP[mem_mat[, g]]),
            all_groups
        ))

        if (length(clump_sets) >= 2) {
            clump_ud <- UpSetR::fromList(clump_sets)

            save_upset(
                suppressWarnings(UpSetR::upset(
                    clump_ud,
                    sets            = rev(names(clump_sets)),
                    order.by        = "freq",
                    mb.ratio        = c(0.65, 0.35),
                    text.scale      = c(1.2, 1.1, 1, 1, 1.1, 0.9),
                    point.size      = 2.8,
                    line.size       = 0.8,
                    mainbar.y.label = sprintf("Overlapping QTLs (±%g kb)", args$overlap_kb),
                    sets.x.label    = "QTLs per group"
                )),
                path = file.path(args$outdir,
                                 paste0("qtl_overlap_upset_clump_", args$name, ".pdf"))
            )
            message("Written: qtl_overlap_upset_clump_", args$name, ".pdf")

            if (length(clump_sets) <= 4) {
                venn_clump <- ggVennDiagram::ggVennDiagram(clump_sets, label_alpha = 0) +
                    scale_fill_gradient(low = "white", high = SIG_STRICT_COLOR) +
                    scale_color_manual(values = rep("grey40", length(clump_sets))) +
                    theme(legend.position = "none")
                save_plot(venn_clump,
                          path   = file.path(args$outdir,
                                             paste0("qtl_overlap_venn_clump_", args$name, ".pdf")),
                          width  = FIG$half_width * 1.5,
                          height = FIG$panel_height * 1.2)
                message("Written: qtl_overlap_venn_clump_", args$name, ".pdf")
            }
        } else {
            message("Fewer than 2 groups with CLUMP QTLs — skipping E2 overlap plots")
        }
    } else {
        message("all.clumped is empty — skipping E2 overlap plots")
    }
}
