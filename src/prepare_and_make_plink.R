#!/usr/bin/env Rscript
#
# prepare_and_make_plink.R
#
# Reads phenotype CSV + YAML config, processes genotypes, filters samples,
# scales phenotypes, and writes GEMMA-format files plus PLINK BED/BIM/FAM.
#
# Outputs (all written to --outdir):
#   pheno_{name}.txt          — phenotype matrix for GEMMA (no header)
#   anno_{name}.txt           — SNP annotations for GEMMA (chr, bp, rs)
#   covars_{name}.txt         — covariate matrix for GEMMA (no header)
#   annotations.csv           — SNP annotations with header (for downstream R)
#   phenotypes_order.txt      — ordered phenotype names (one per line)
#   raw_phenotypes.csv        — unscaled phenotypes with Strain/MouseID
#   export_strains_order.csv  — ordered strain names after downsampling
#   {name}.bed/.bim/.fam      — PLINK files for GEMMA kinship + SuSiE
#
# Usage:
#   Rscript prepare_and_make_plink.R \
#     --input      data.csv         \
#     --yaml       run.yaml         \
#     --name       my_run           \
#     --genotypes  hmdp.csv.gz      \
#     --outdir     results/input

library(argparse)
library(yaml)
library(dplyr)
library(readr)
library(data.table)
library(genio)

# ── average_strain ────────────────────────────────────────────────────────────
# Inlined from mousegwas/R/gemma_tools.R to remove package dependency.
# Returns a list with $genotypes, $phenotypes, $indices, $covars after
# downsampling or averaging within each strain.
average_strain <- function(strains_genomes, phenotypes, covars, downsample, sex, strain = NULL) {
    set.seed(100)
    cret <- NULL
    if (is.null(strain)) {
        sgs   <- rbind(as.list(c(1, 1, 1, sex)), strains_genomes, use.names = FALSE)
        grows <- c(1, sample(nrow(sgs), min(nrow(sgs), 1000)))
        genidx <- match(sgs[grows, ], sgs[grows, ])
    } else {
        sgs    <- paste(c(1, 2, 3, strain), c(1, 2, 3, sex))
        genidx <- match(sgs, sgs)
    }

    if (downsample > 0) {
        miceidx <- 1:3   # 1:3 = rs, minor, major columns
        for (i in unique(genidx[-1:-3])) {
            idx <- which(genidx == i)
            miceidx <- c(miceidx, if (length(idx) > 1) sample(idx, min(downsample, length(idx))) else idx)
        }
    } else {
        miceidx <- which(!duplicated(genidx))
    }

    gret <- strains_genomes[, miceidx, with = FALSE]

    if (downsample == 0) {
        phen2         <- cbind(as.data.frame(phenotypes), as.data.frame(covars[, -1, drop = FALSE]))
        phen2$strain  <- factor(genidx[-1:-3])
        pret          <- NULL
        for (pn in colnames(phenotypes)) {
            lmout <- lm(as.formula(paste0(
                pn, " ~ 0 + ", paste(colnames(covars)[-1], collapse = "+"), " + strain"
            )), phen2)
            pret <- cbind(pret, lmout$coefficients[grepl("^strain", names(lmout$coefficients))])
        }
        pret <- setNames(as.data.frame(pret), colnames(phenotypes))
    } else {
        pret <- phenotypes[miceidx[-1:-3] - 3, , drop = FALSE]
        cret <- if (!is.null(covars)) covars[miceidx[-1:-3] - 3, , drop = FALSE] else NULL
    }

    list(genotypes = gret, phenotypes = pret, indices = miceidx[-1:-3] - 3, covars = cret)
}

# ── Argument parsing ──────────────────────────────────────────────────────────
parser <- ArgumentParser()

parser$add_argument("--input",      "-i", required = TRUE,
    help = "Input phenotype CSV (column names required)")
parser$add_argument("--yaml",       "-y", required = TRUE,
    help = "YAML config: strain, phenotypes, covariates, F1 mappings, etc.")
parser$add_argument("--name",       "-n", required = TRUE,
    help = "Run name — used as prefix for all output files")
parser$add_argument("--genotypes",  "-g", required = TRUE,
    help = "HMDP genotype CSV (e.g. muster_hmdp.csv.gz)")
parser$add_argument("--outdir",     "-o", default = ".",
    help = "Output directory (created if absent, default: .)")
parser$add_argument("--downsample", "-d", default = 100000L, type = "integer",
    help = "Max individuals per strain; 0 = average within strain (default: 100000)")
parser$add_argument("--MAF",     default = 0.1,  type = "double",
    help = "Min minor allele frequency (default: 0.1)")
parser$add_argument("--missing", default = 0.05, type = "double",
    help = "Max fraction of missing genotype calls per SNP (default: 0.05)")
parser$add_argument("--qqnorm",  default = FALSE, action = "store_true",
    help = "QQ-normalise each phenotype before output")
parser$add_argument("--shuffle", default = FALSE, action = "store_true",
    help = "Permute phenotype labels (for null-distribution runs)")
parser$add_argument("--seed",    default = 100L, type = "integer",
    help = "RNG seed for --shuffle (default: 100)")

args <- parser$parse_args()

# ── Setup ─────────────────────────────────────────────────────────────────────
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
name <- args$name

# ── Load YAML and phenotype table ─────────────────────────────────────────────
yamin          <- yaml.load_file(args$yaml)
complete_table <- read_csv(args$input, col_types = cols(), show_col_types = FALSE)

# ── Dataset-specific filters ──────────────────────────────────────────────────
# NOTE: these filters are specific to the anxfb coat-colour analysis.
complete_table <- complete_table |> mutate(CoatColorGWAS = CoatColor) 

# ── Build phenotype and covariate name lists from YAML ────────────────────────
pheno_names <- c()
phegroups   <- list()
i <- 1L
for (n in names(yamin$phenotypes)) {
    if (n %in% names(complete_table)) {
        pheno_names <- c(pheno_names, n)
        grp <- yamin$phenotypes[[n]]$group
        if (!is.null(grp))
            phegroups[[grp]] <- c(phegroups[[grp]], i)
        i <- i + 1L
    }
}

covar_names <- setdiff(intersect(yamin$covar, names(complete_table)), pheno_names)

# ── Resolve strain names (inbreds + F1s + translations) ──────────────────────
strains <- complete_table |>
    select(input_name = one_of(yamin$strain)) |>
    distinct() |>
    mutate(p1 = input_name, p2 = input_name)

f1_map <- tibble(input_name = character(), p1 = character(), p2 = character())
for (fl in yamin$F1)
    f1_map <- add_row(f1_map, input_name = names(fl), p1 = unlist(fl)[1], p2 = unlist(fl)[2])

strains <- strains |>
    left_join(f1_map, by = "input_name", suffix = c(".orig", "")) |>
    mutate(p1 = coalesce(p1, p1.orig), p2 = coalesce(p2, p2.orig)) |>
    select(input_name, p1, p2)

trans <- tibble(name = character(), gen_name = character())
for (tr in yamin$translate)
    trans <- add_row(trans, name = names(tr), gen_name = unlist(tr)[1])

strains <- strains |>
    left_join(trans, by = c("p1" = "name")) |>
    mutate(p1 = coalesce(gen_name, p1)) |>
    select(input_name, p1, p2) |>
    left_join(trans, by = c("p2" = "name")) |>
    mutate(p2 = coalesce(gen_name, p2)) |>
    select(input_name, p1, p2)

valid_strains <- unique(c(strains$p1, strains$p2))

# ── Load and recode genotypes ─────────────────────────────────────────────────
message("Loading genotypes from: ", args$genotypes)
geno <- fread(args$genotypes)
geno[, c("major", "minor") := tstrsplit(observed, "/", fixed = TRUE, keep = 1:2)]
geno <- geno[rs != ""]

keep_cols   <- c(intersect(names(geno), valid_strains), "chr", "bp38", "rs", "major", "minor")
complete.geno <- geno[, ..keep_cols]
complete.geno[, chr := as.character(chr)]

for (cn in setdiff(names(complete.geno), c("chr", "bp38", "rs", "major", "minor"))) {
    complete.geno[get(cn) == "H",   (cn) := 1]
    complete.geno[get(cn) == major, (cn) := 0]
    complete.geno[get(cn) == minor, (cn) := 2]
    complete.geno[, (cn) := as.numeric(get(cn))]
}

srdata <- complete.geno[, .(rs, major, minor)]

# ── Build per-individual genotype, phenotype, and covariate tables ────────────
strains_genomes <- srdata
sorder   <- c()
phenos   <- data.table()
for (p in pheno_names) phenos[, (p) := numeric()]
covars   <- data.table()
for (cv in covar_names) covars[, (cv) := numeric()]
covars[, isWild := numeric()]
sexvec   <- c()
notfound <- c()

for (comrow in seq_len(nrow(complete_table))) {
    sname <- as.character(complete_table[[comrow, yamin$strain]])
    rnum  <- which(strains$input_name == sname)
    if (length(rnum) == 0) {
        message("Strain not found in config: ", sname, " (row ", comrow, ")")
        next
    }
    p1n <- strains$p1[rnum]
    p2n <- strains$p2[rnum]

    if (!(p1n %in% names(complete.geno)) || !(p2n %in% names(complete.geno))) {
        key <- paste(p1n, p2n)
        if (!key %in% notfound) {
            message("Genotype not found: ", p1n, " / ", p2n)
            notfound <- c(notfound, key)
        }
        next
    }

    sorder <- c(sorder, sname)
    strains_genomes[, paste0("X", comrow) :=
        (complete.geno[[p1n]] + complete.geno[[p2n]]) / 2]

    # X chromosome in F1 males: take only the maternal haplotype
    if (p1n != p2n && !is.null(yamin$sex) && complete_table[[comrow, yamin$sex]] == "M")
        strains_genomes[complete.geno$chr == "X", paste0("X", comrow) :=
            complete.geno[[p1n]][complete.geno$chr == "X"]]

    prow <- complete_table[comrow, c("Strain", "MouseID", pheno_names)]
    crow <- cbind(
        complete_table[comrow, covar_names, drop = FALSE],
        tibble(isWild = as.numeric(p1n %in% yamin$wild | p2n %in% yamin$wild))
    )

    if (length(covar_names) == 0 || all(!is.na(crow))) {
        phenos <- rbind(phenos, prow, fill = TRUE)
        covars <- rbind(covars, crow, fill = TRUE)
        sexvec <- c(sexvec,
            if (is.null(yamin$sex)) 1L else complete_table[[comrow, yamin$sex]])
    }
}

# ── Build covariate model matrix ──────────────────────────────────────────────
if (length(covar_names) > 0) {
    covars <- model.matrix(
        as.formula(paste0("~", paste(covar_names, collapse = "+"))),
        data = covars
    )
} else {
    covars <- NULL
}

# ── Scale continuous phenotypes ───────────────────────────────────────────────
raw_phenos <- copy(phenos)

cols_to_scale <- sapply(phenos, function(x)
    !is.factor(x) && !is.character(x) && !all(is.na(x) | x == 0 | x == 1))
scale_cols <- names(which(cols_to_scale))
phenos[, (scale_cols) := lapply(.SD, scale), .SDcols = scale_cols]

# Drop all-NA or all-zero columns
for (cc in rev(seq_len(ncol(phenos)))) {
    if (all(is.na(phenos[[cc]]) | phenos[[cc]] == 0)) {
        phenos     <- phenos[, -cc, with = FALSE]
        raw_phenos <- raw_phenos[, -cc, with = FALSE]
    }
}

# Optional permutation
if (args$shuffle) {
    set.seed(args$seed)
    neword     <- sample(nrow(phenos))
    phenos     <- phenos[neword, ]
    raw_phenos <- raw_phenos[neword, ]
}

# ── Strain averaging / downsampling ──────────────────────────────────────────
b          <- average_strain(strains_genomes, phenos, covars, args$downsample, sexvec, sorder)
raw_phenos <- raw_phenos[b$indices, , drop = FALSE]

# ── QQ-normalise (optional) ───────────────────────────────────────────────────
if (args$qqnorm)
    for (r in seq_len(ncol(b$phenotypes)))
        b$phenotypes[, r] <- qqnorm(as.data.frame(b$phenotypes)[, r], plot = FALSE)$x

# ── SNP filters: MAF, missingness, heterozygous calls ─────────────────────────
n_ind <- ncol(complete.geno) - 5L
mafc  <- rowSums(complete.geno[, !c("chr", "bp38", "rs", "major", "minor")],
                 na.rm = TRUE) / (2 * n_ind)

b$genotypes <- b$genotypes[
    rowSums(is.na(b$genotypes)) <= (ncol(b$genotypes) - 3L) * args$missing &
    mafc >= args$MAF & mafc <= 1 - args$MAF
]

het_rs      <- b$genotypes[
    rowSums(b$genotypes[, !c("rs", "major", "minor")] == 0.5, na.rm = TRUE) > 0, rs]
b$genotypes <- b$genotypes[!rs %in% het_rs]

complete.geno <- complete.geno[rs %in% b$genotypes$rs]

# ── Write GEMMA-format outputs ────────────────────────────────────────────────
anno_cols   <- c("chr", "bp38", "rs")
annotations <- complete.geno[, ..anno_cols]

fwrite(annotations,
    file.path(args$outdir, paste0("anno_", name, ".txt")),
    col.names = FALSE, sep = "\t")

fwrite(annotations,
    file.path(args$outdir, "annotations.csv"),
    col.names = TRUE, sep = ",")

pheno_mat <- b$phenotypes[, !colnames(b$phenotypes) %in% c("Strain", "MouseID"), drop = FALSE]
fwrite(as.data.table(pheno_mat),
    file.path(args$outdir, paste0("pheno_", name, ".txt")),
    col.names = FALSE, sep = "\t")

if (!is.null(b$covars))
    fwrite(as.data.table(b$covars),
        file.path(args$outdir, paste0("covars_", name, ".txt")),
        col.names = FALSE, sep = "\t")

write.csv(
    colnames(b$phenotypes)[!colnames(b$phenotypes) %in% c("Strain", "MouseID")],
    file.path(args$outdir, "phenotypes_order.txt"),
    quote = FALSE, col.names = FALSE, row.names = FALSE)

write.csv(raw_phenos,
    file.path(args$outdir, "raw_phenotypes.csv"),
    row.names = FALSE)

write.table(sorder[b$indices],
    file.path(args$outdir, "export_strains_order.csv"),
    quote = FALSE, col.names = FALSE)

# ── Write PLINK BED/BIM/FAM ───────────────────────────────────────────────────
message("Writing PLINK files...")

plink_pheno     <- as.data.frame(b$phenotypes)
plink_pheno$FID <- as.integer(as.factor(plink_pheno$Strain))
plink_pheno$IID <- seq_len(nrow(plink_pheno))

m <- nrow(b$genotypes)
n <- ncol(b$genotypes) - 3L   # exclude rs, major, minor

X <- as.matrix(b$genotypes[, !c("rs", "major", "minor")])
rownames(X) <- b$genotypes$rs
colnames(X) <- seq_len(n)

bim      <- make_bim(n = m)
bim$chr  <- annotations$chr
bim$id   <- annotations$rs
bim$pos  <- annotations$bp38
bim$ref  <- ifelse(is.na(b$genotypes$minor) | b$genotypes$minor == "", "0", b$genotypes$minor)
bim$alt  <- ifelse(is.na(b$genotypes$major) | b$genotypes$major == "", "0", b$genotypes$major)

fam        <- make_fam(n = n)
fam$fam    <- plink_pheno$FID
fam$id     <- plink_pheno$IID
fam$sex    <- ifelse(complete_table$Sex[b$indices] == "M", 1L, 2L)
fam$pheno  <- -9L

# Drop sex chromosomes (avoids hemizygosity warnings in PLINK/GEMMA)
autosomal  <- !bim$chr %in% c("X", "Y", "MT")
bim        <- bim[autosomal, ]
X          <- X[autosomal, ]

# Sort by chromosome then position (required by PLINK)
sort_idx   <- order(as.integer(bim$chr), bim$pos)
bim        <- bim[sort_idx, ]
X          <- X[sort_idx, ]

write_plink(file.path(args$outdir, name), X, bim, fam)

message("Done. Outputs written to: ", normpath(args$outdir))
