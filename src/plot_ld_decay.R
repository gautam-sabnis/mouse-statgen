#!/usr/bin/env Rscript
#
# plot_ld_decay.R
#
# Reads pairwise r² from PLINK --r2 output, bins pairs by physical distance,
# and plots genome-wide LD decay with loess smoothing.
#
# Input:
#   <prefix>.ld  — PLINK --r2 output (CHR_A BP_A SNP_A CHR_B BP_B SNP_B R^2)
#
# Output:
#   LD_decay_<name>.pdf

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

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
})

if (LOCAL) {
    args <- list(
        ld      = "output/input/ld_all.ld",
        name    = "anxfb_test",
        outdir  = "output/plots",
        binsize = 5000L,
        maxdist = 2.5
    )
} else {
    parser <- ArgumentParser()
    parser$add_argument("--ld", required = TRUE,
        help = "PLINK .ld file from --r2 computation")
    parser$add_argument("--name", "-n", required = TRUE,
        help = "Run name — used in output filename")
    parser$add_argument("--outdir", "-o", default = ".",
        help = "Output directory (default: current dir)")
    parser$add_argument("--binsize", type = "integer", default = 5000L,
        help = "Bin width in bp for averaging r² (default: 5000)")
    parser$add_argument("--maxdist", type = "double", default = 2.5,
        help = "Maximum distance to plot in Mbp (default: 2.5)")
    args <- parser$parse_args()
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# Read PLINK .ld — select R^2 by position (last column) to avoid caret parsing issues
ld_raw <- fread(args$ld)
r2_col <- tail(names(ld_raw), 1L)   # always the last column in PLINK .ld output
ld <- data.table(
    dist = abs(ld_raw$BP_A - ld_raw$BP_B),
    r_sq = ld_raw[[r2_col]]
)
rm(ld_raw)

# Bin by distance, average r² within each bin
ld[, distc := cut(dist,
    breaks = seq(0, max(dist) + args$binsize, by = args$binsize),
    include.lowest = TRUE)]
avg <- ld[!is.na(distc),
          .(avdist = mean(dist), avr_sq = mean(r_sq, na.rm = TRUE)),
          by = distc]

plt <- ggplot(avg, aes(avdist / 1e6, avr_sq)) +
    geom_smooth(
        color  = LD_COLOR,
        method = "loess",
        span   = 0.3,
        se     = FALSE
    ) +
    xlim(c(0, args$maxdist)) +
    labs(x = "Distance (Mbp)",
         y = expression("Average LD" ~ (r^2))) +
    theme_gwas()

out <- file.path(args$outdir, paste0("LD_decay_", args$name, ".pdf"))
save_plot(plt, path = out, width = FIG$half_width, height = FIG$panel_height)
message("Written: ", out)
