#!/usr/bin/env Rscript
#
# aesthetics.R
#
# Central aesthetic configuration for all postprocessing scripts.
# Source this file once to get consistent colors, themes, and dimensions.
# To restyle all figures at once, change values here — nowhere else.
#
# Requires: ggplot2, ggsci, viridis
# ggsci provides publication-quality categorical palettes (NPG, NEJM, BMJ, JCO,
# Lancet, JAMA) as both raw color vectors and ready-to-use ggplot2 scales.

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggsci)
    library(viridis)
})

# ── Figure dimensions (inches, publication layout) ────────────────────────────

FIG <- list(
    full_width   = 7.25,
    half_width   = 3.54,
    full_height  = 9.75,    # full page minus margins (11 - 1.25)
    panel_height = 3.54     # height of a single panel row
)

# Font: falls back to the device default ("") if Arial is not registered locally.
# On HPC, register Arial first with extrafont::loadfonts() or equivalent.
FIG_FONT <- tryCatch({ grDevices::pdfFonts("Arial"); "Arial" },
                     error = function(e) { message("Arial not found, using default font"); "" })

# ── Color palettes ─────────────────────────────────────────────────────────────
# All categorical palettes come from ggsci. The palette argument accepts any of:
#   "npg", "nejm", "bmj", "jco", "lancet", "jama", "d3"
# Viridis is kept for continuous/sequential scales (heatmaps etc.).

# Cluster colors: discrete, length determined at call time by k.
# Default IGV — 51 colors from the Integrative Genomics Viewer palette,
# designed for genomics and maximally distinct even at high k.
# D3 (10 colors) is a good alternative for k <= 10.
cluster_colors <- function(k, palette = c("igv", "d3", "npg", "jco", "nejm", "lancet")) {
    palette <- match.arg(palette)
    pal_fn  <- switch(palette,
        igv    = ggsci::pal_igv(),
        d3     = ggsci::pal_d3(),
        npg    = ggsci::pal_npg(),
        jco    = ggsci::pal_jco(),
        nejm   = ggsci::pal_nejm(),
        lancet = ggsci::pal_lancet()
    )
    pal_fn(k)
}

# Group colors: named vector, one color per group label.
# Default NEJM (8 colors); swap palette as needed.
group_colors <- function(groups,
                         palette = c("nejm", "npg", "bmj", "jco", "lancet", "jama")) {
    palette <- match.arg(palette)
    pal_fn  <- switch(palette,
        nejm   = ggsci::pal_nejm(),
        npg    = ggsci::pal_npg(),
        bmj    = ggsci::pal_bmj(),
        jco    = ggsci::pal_jco(),
        lancet = ggsci::pal_lancet(),
        jama   = ggsci::pal_jama()
    )
    stats::setNames(pal_fn(length(groups)), groups)
}

# Heatmap color scale (peaks x phenotypes matrix) — continuous, keep viridis
HEATMAP_COLORS <- viridis::viridis(128)

# Chromosome background: alternating grey for Manhattan plots
CHROM_COLORS <- rep(c("#CCCCCC", "#969696"), 10)

# Significance threshold reference line
SIG_LINE_COLOR <- "#FCBBA1"

# Significant SNP point colors (Manhattan plots)
SIG_STRICT_COLOR  <- "#CC2200"   # above the stricter threshold
SIG_LENIENT_COLOR <- "#E69F00"   # above the lenient threshold only

# MAF histogram: [1] all SNPs, [2] peak SNPs — first two NEJM colors
MAF_COLORS <- ggsci::pal_nejm()(2)

# LD decay smoothed line — third JCO color
LD_COLOR <- ggsci::pal_jco()(3)[3]

# ── ggsci scale shortcuts ──────────────────────────────────────────────────────
# Pre-built ggplot2 scales; add directly to a plot with `+`.
# Example:  p + SCALE_COLOR_NPG   or   p + SCALE_FILL_NEJM

SCALE_COLOR_NPG    <- ggsci::scale_color_npg()
SCALE_FILL_NPG     <- ggsci::scale_fill_npg()
SCALE_COLOR_NEJM   <- ggsci::scale_color_nejm()
SCALE_FILL_NEJM    <- ggsci::scale_fill_nejm()
SCALE_COLOR_BMJ    <- ggsci::scale_color_bmj()
SCALE_FILL_BMJ     <- ggsci::scale_fill_bmj()
SCALE_COLOR_JCO    <- ggsci::scale_color_jco()
SCALE_FILL_JCO     <- ggsci::scale_fill_jco()
SCALE_COLOR_LANCET <- ggsci::scale_color_lancet()
SCALE_FILL_LANCET  <- ggsci::scale_fill_lancet()
SCALE_COLOR_JAMA   <- ggsci::scale_color_jama()
SCALE_FILL_JAMA    <- ggsci::scale_fill_jama()

# ── ggplot2 themes ─────────────────────────────────────────────────────────────

# Base theme: bw with no border, no x-axis gridlines
theme_gwas <- function(base_size = 10, family = FIG_FONT) {
    theme_bw(base_size = base_size) %+replace%
        theme(
            text               = element_text(size = base_size, family = family),
            panel.border       = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank()
        )
}

# Manhattan theme: extends base, drops legend and y gridlines
theme_manhattan <- function(base_size = 10, family = FIG_FONT) {
    theme_gwas(base_size, family) %+replace%
        theme(
            legend.position    = "none",
            panel.grid.major.y = element_blank()
        )
}

# PVE bar chart theme: extends base, legend on right
theme_pve <- function(base_size = 10, family = FIG_FONT) {
    theme_gwas(base_size, family) %+replace%
        theme(
            axis.ticks.y    = element_blank(),
            legend.position = "right"
        )
}

# General-purpose minimal theme: clean, base-14, for reports and exploratory plots
theme_base <- function(base_size = 14, family = FIG_FONT) {
    theme_minimal(base_size = base_size) %+replace%
        theme(
            text         = element_text(size = base_size, family = family),
            axis.text.x  = element_text(size = base_size),
            axis.text.y  = element_text(size = base_size),
            legend.title = element_text(size = base_size),
            legend.text  = element_text(size = base_size)
        )
}

# Pre-built shorthand — use `mytheme` directly in any postprocessing script
mytheme <- theme_base()

# ── Plot saving wrapper ────────────────────────────────────────────────────────
# Thin wrapper around ggsave with publication defaults.
# Prefers cairo_pdf for font rendering; falls back to pdf if Cairo is not
# available (common on local macOS installs without a Cairo-enabled R build).
# Pass width / height explicitly to override FIG defaults.

PLOT_DEVICE <- if (capabilities("cairo")) cairo_pdf else {
    message("Cairo not available, falling back to pdf device")
    grDevices::pdf
}

save_plot <- function(plot, path,
                      width  = FIG$full_width,
                      height = FIG$panel_height,
                      device = PLOT_DEVICE,
                      ...) {
    ggplot2::ggsave(
        filename = path,
        plot     = plot,
        device   = device,
        dpi      = "print",
        width    = width,
        height   = height,
        units    = "in",
        ...
    )
}
