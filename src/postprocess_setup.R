#!/usr/bin/env Rscript
#
# postprocess_setup.R
#
# Data loading functions shared across all postprocessing scripts.
# Source this file at the top of each analysis script:
#
#   source(file.path(SRC, "aesthetics.R"))       # must come first
#   source(file.path(SRC, "postprocess_setup.R"))
#
# where SRC is the directory containing both files (typically projectDir/src).
#
# Functions:
#   load_phenotype_config()  -- phenotype names, groups, paper names from YAML
#   assign_group_colors()    -- attach group color column to pnames
#   load_annotations()       -- SNP chromosome/position table
#   load_genotypes()         -- genotype matrix (tidy tibble + numeric matrix)
#   load_pve()               -- PVE estimates joined with phenotype metadata

if (!exists("theme_gwas"))
    stop("Source aesthetics.R before postprocess_setup.R")

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tibble)
    library(yaml)
})

# ── Phenotype configuration ────────────────────────────────────────────────────
# Parses the YAML file and builds a phenotype metadata table.
#
# Arguments:
#   yaml_file      -- path to the pipeline YAML
#   coat_phenotype -- if TRUE, add coat/eye phenotypes from phenos_order
#   phenos_order   -- character vector of phenotype IDs (from phenotypes_order.txt);
#                     required when coat_phenotype = TRUE
#
# Returns a list:
#   $pnames       -- data.frame(PaperName, Group), rownames = internal phenotype names
#   $pheno_names  -- character vector of internal phenotype names (YAML order)
#   $groups_order -- character vector of group names
#                    (from yaml$groups key, or unique(pnames$Group) if absent)
#   $yaml         -- the raw parsed YAML object (for downstream use)
load_phenotype_config <- function(yaml_file,
                                  coat_phenotype = FALSE,
                                  phenos_order   = NULL) {
    yamin <- yaml::yaml.load_file(yaml_file)

    pnames      <- data.frame(PaperName = character(0),
                              Group     = character(0),
                              stringsAsFactors = FALSE)
    pheno_names <- character(0)

    for (n in names(yamin$phenotypes)) {
        ph          <- yamin$phenotypes[[n]]
        pheno_names <- c(pheno_names, n)
        pnames      <- rbind(pnames, data.frame(
            Group     = if (length(ph$group)) ph$group else "NoGroup",
            PaperName = ph$papername,
            stringsAsFactors = FALSE
        ))
    }

    if (coat_phenotype && !is.null(phenos_order)) {
        for (p in phenos_order) {
            pheno_names <- c(pheno_names, p)
            pnames      <- rbind(pnames, data.frame(
                Group     = if (grepl("coat", p)) "Coat" else "Eyes",
                PaperName = gsub(" ", "-", gsub("(coat|eyes)", "", p)),
                stringsAsFactors = FALSE
            ))
        }
    }

    rownames(pnames) <- pheno_names

    groups_order <- if (length(yamin$groups)) yamin$groups
                    else                       unique(pnames$Group)

    list(pnames       = pnames,
         pheno_names  = pheno_names,
         groups_order = groups_order,
         yaml         = yamin)
}

# ── Group color assignment ─────────────────────────────────────────────────────
# Attaches a 'color' column to pnames based on Group membership.
# Relies on group_colors() from aesthetics.R.
#
# meanvariance mode: variance groups inherit the color of their mean counterpart
# (e.g. "WeightVariance" gets the same color as "Weight").
#
# Arguments:
#   pnames       -- data.frame from load_phenotype_config()$pnames
#   groups_order -- character vector from load_phenotype_config()$groups_order
#   palette      -- ggsci palette name passed to group_colors()
#                   ("nejm", "npg", "bmj", "jco", "lancet", "jama")
#   meanvariance -- logical; pair variance groups with their mean group color
#
# Returns a list:
#   $pnames  -- pnames with 'color' column added
#   $grpcol  -- named character vector (group name -> hex color)
assign_group_colors <- function(pnames,
                                groups_order,
                                palette      = "nejm",
                                meanvariance = FALSE) {
    if (meanvariance) {
        mean_groups <- groups_order[!grepl("Variance", groups_order)]
        grpcol      <- group_colors(mean_groups, palette = palette)
        for (g in groups_order[grepl("Variance", groups_order)]) {
            base      <- gsub("Variance", "", g)
            grpcol[g] <- unname(grpcol[base])
        }
    } else {
        grpcol <- group_colors(groups_order, palette = palette)
    }

    pnames <- dplyr::left_join(
        pnames,
        tibble::tibble(Group = names(grpcol), color = unname(grpcol)),
        by = "Group"
    )
    rownames(pnames) <- rownames(pnames)  # left_join drops rownames; restore below
    # rownames were lost by left_join — caller must restore if needed
    # (use: rownames(result$pnames) <- cfg$pheno_names)

    list(pnames = pnames, grpcol = grpcol)
}

# ── SNP annotations ────────────────────────────────────────────────────────────
# Reads annotations.csv produced by the pipeline.
#
# Returns a tibble with columns: rs (character), ps (double), chr (character)
load_annotations <- function(outdir) {
    readr::read_delim(
        file.path(outdir, "annotations.csv"),
        delim          = ",",
        col_names      = c("rs", "ps", "chr"),
        col_types      = readr::cols(rs  = readr::col_character(),
                                     ps  = readr::col_double(),
                                     chr = readr::col_character()),
        show_col_types = FALSE
    )
}

# ── Genotype matrix ────────────────────────────────────────────────────────────
# Reads strains_genotypes_all.csv.
# The tidy form (geno_t) is used for MAF / LD calculations that need metadata.
# The matrix form (geno) is passed to mousegwas peak-calling functions.
#
# Returns a list:
#   $geno_t -- tibble: rs, chr, bp38, major, minor, <one column per strain>
#   $geno   -- numeric matrix: rownames = rs IDs, columns = strain IDs
load_genotypes <- function(outdir) {
    geno_t <- readr::read_csv(
        file.path(outdir, "strains_genotypes_all.csv"),
        col_types      = readr::cols(.default = readr::col_double(),
                                     chr      = readr::col_character(),
                                     rs       = readr::col_character(),
                                     major    = readr::col_character(),
                                     minor    = readr::col_character()),
        show_col_types = FALSE
    )
    geno <- as.matrix(
        geno_t |>
            tibble::column_to_rownames("rs") |>
            dplyr::select(-chr, -bp38, -major, -minor)
    )
    list(geno_t = geno_t, geno = geno)
}

# ── PVE estimates ──────────────────────────────────────────────────────────────
# Reads PVE_GEMMA_estimates.txt and left-joins Group / PaperName from pnames.
# pnames must have rownames = internal phenotype names.
#
# Returns a tibble with all PVE columns plus Group, PaperName, color.
load_pve <- function(outdir, pnames) {
    pve <- readr::read_csv(
        file.path(outdir, "PVE_GEMMA_estimates.txt"),
        show_col_types = FALSE
    )
    dplyr::left_join(
        pve,
        tibble::as_tibble(pnames, rownames = "phenotype"),
        by = "phenotype"
    )
}
