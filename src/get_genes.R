#!/usr/bin/env Rscript
#
# get_genes.R
#
# Annotates GWAS results with nearby mouse genes via Ensembl biomaRt (GRCm38).
#
# Input modes (--mode):
#   susie  — annotate credible set SNPs from SuSiE output files (PIPs + CS)
#   clump  — annotate lead SNPs from PLINK clump file
#
# Outputs:
#   gene_annotations.csv  — SNP-gene pairs with window overlaps

library(argparse)
library(data.table)
library(biomaRt)

parser <- ArgumentParser()
parser$add_argument("--mode",    required = TRUE,
    help = "Input mode: 'susie' (credible set SNPs) or 'clump' (lead SNPs)")
parser$add_argument("--input",   required = TRUE,  nargs = "+",
    help = "Input file(s): SuSiE locus txt files or PLINK .clumped file")
parser$add_argument("--dist",    default  = 1e6,   type = "double",
    help = "Window in bp around each SNP to search for genes (default 1Mb)")
parser$add_argument("--host",    default  = "https://sep2025.archive.ensembl.org",
    help = "Ensembl archive host (default: GRCm38 / Ensembl 102)")
parser$add_argument("--attempts",default  = 5L,    type = "integer",
    help = "Number of biomaRt retry attempts (default 5)")
parser$add_argument("--outdir",  default  = ".",
    help = "Output directory (default: current dir)")
args <- parser$parse_args()

stopifnot(args$mode %in% c("susie", "clump"))

# ── Build SNP table from input ────────────────────────────────────────────────
if (args$mode == "susie") {
    # Collect all credible set SNPs (cs not NA) across all locus files
    snps <- rbindlist(lapply(args$input, function(f) {
        dt <- fread(f)
        dt[!is.na(cs), .(chr = as.character(chr), ps = bp, snp, pip, cs, trait, lead_snp)]
    }), fill = TRUE)
    if (nrow(snps) == 0) stop("No credible set SNPs found in SuSiE input files.")

} else {
    # Use lead SNPs from PLINK clump file
    snps <- rbindlist(lapply(args$input, function(f) {
        dt <- fread(f)
        if (nrow(dt) == 0) return(NULL)
        dt[, .(chr = as.character(CHR), ps = BP, snp = SNP)]
    }), fill = TRUE)
    if (nrow(snps) == 0) stop("No lead SNPs found in clump input files.")
}

cat(sprintf("Mode: %s | %d SNPs to annotate\n", args$mode, nrow(snps)))

# ── Fetch gene annotations from Ensembl ───────────────────────────────────────
annot  <- NULL
attempt <- 0L
while (attempt < args$attempts && is.null(annot)) {
    attempt <- attempt + 1L
    ensembl <- NULL
    try(silent = TRUE,
        ensembl <- biomaRt::useMart("ensembl",
                                    dataset = "mmusculus_gene_ensembl",
                                    host    = args$host))
    if (is.null(ensembl)) next
    try(silent = TRUE,
        annot <- biomaRt::getBM(
            attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name",
                           "strand", "start_position", "end_position", "gene_biotype"),
            mart     = ensembl,
            useCache = TRUE
        ))
}
if (is.null(annot))
    stop(sprintf("biomaRt query failed after %d attempts.", args$attempts))

cat(sprintf("Retrieved %d genes from Ensembl.\n", nrow(annot)))

# ── Join SNPs to genes within window ─────────────────────────────────────────
annot_dt <- as.data.table(annot)
annot_dt <- annot_dt[!is.na(chromosome_name) &
                     !is.na(start_position)   &
                     !is.na(end_position)]

# Pre-filter genes to only chromosomes present in the SNP table
keep_chrs <- unique(snps$chr)
annot_dt  <- annot_dt[chromosome_name %in% keep_chrs]
annot_dt[, ext_start := start_position - args$dist]
annot_dt[, ext_end   := end_position   + args$dist]

# Join per chromosome to avoid memory explosion from full cross-join
result <- rbindlist(lapply(keep_chrs, function(ch) {
    snps_ch  <- snps[chr == ch]
    annot_ch <- annot_dt[chromosome_name == ch]
    if (nrow(snps_ch) == 0 || nrow(annot_ch) == 0) return(NULL)
    # Range join: SNP position falls within extended gene window
    snps_ch[, ps_end := ps]   # point SNP — start == end
    setkey(annot_ch, ext_start, ext_end)
    hits <- foverlaps(snps_ch, annot_ch,
                      by.x    = c("ps", "ps_end"),
                      by.y    = c("ext_start", "ext_end"),
                      type    = "within",
                      nomatch = NA)
    hits[, ps_end := NULL]
    hits
}), fill = TRUE)

cat(sprintf("Annotated %d SNP-gene pairs.\n", nrow(result)))

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
fwrite(result, file.path(args$outdir, "gene_annotations.csv"))
cat(sprintf("Written to %s/gene_annotations.csv\n", args$outdir))
