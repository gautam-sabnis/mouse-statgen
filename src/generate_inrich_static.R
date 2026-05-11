#!/usr/bin/env Rscript
#
# generate_inrich_static.R
#
# One-time script to pre-generate all static INRICH gene-set input files.
# Run interactively once (requires internet for biomaRt).
# Commit outputs to data/inrich_static/ and point the pipeline at them.
#
# Usage:
#   Rscript src/generate_inrich_static.R \
#       --extdata  data/extdata \
#       --outdir   data/inrich_static
#       --no-sfari  # skip SFARI gene sets (requires human biomaRt access)    
#
# Outputs (all tab-separated, no header unless noted):
#   genes_coordinates_for_INRICH.txt         chr  start  end  ensembl_id  symbol
#   groups_KEGG_pathway_link_for_INRICH.txt  ensembl_id  path_id  path_desc
#   groups_GO_BP_terms_link_for_INRICH.txt   ensembl_id  go_id  go_desc
#   groups_GO_CC_terms_link_for_INRICH.txt   ensembl_id  go_id  go_desc
#   groups_GO_MF_terms_link_for_INRICH.txt   ensembl_id  go_id  go_desc
#   groups_MP_terms_for_INRICH.txt           ensembl_id  mp_id  mp_desc
#   groups_MPMGI_terms_for_INRICH.txt        ensembl_id  mp_id  mp_desc
#   groups_brain_expression_link_for_INRICH.txt  ensembl_id  group_short  group_full
#   groups_SFARI_terms_for_INRICH.txt        ensembl_id  set_id  set_desc

suppressPackageStartupMessages({
    library(argparse)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(tibble)
    library(biomaRt)
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    library(GO.db)
})

parser <- ArgumentParser()
parser$add_argument("--extdata", required = TRUE,
    help = "Path to directory containing the mousegwas extdata files")
parser$add_argument("--outdir", default = "data/inrich_static",
    help = "Output directory for generated files [default: data/inrich_static]")
parser$add_argument("--no-sfari", action = "store_true", default = FALSE,
    help = "Skip SFARI gene-set generation (requires a second biomaRt call for human-mouse ortholog mapping)")
args <- parser$parse_args()

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

write_tsv_nohdr <- function(df, path) {
    write.table(df, file = path, sep = "\t", col.names = FALSE,
                row.names = FALSE, quote = FALSE)
}

# ── 1. biomaRt: fetch all mouse gene annotations (GRCm38) ────────────────────
# Hard-coded to GRCm38 archive — same coordinate system as the MDA genotype chip.
cat("Querying biomaRt (GRCm38)...\n")
ensembl <- NULL
attempts <- 0
while (is.null(ensembl) && attempts < 5) {
    attempts <- attempts + 1
    try(silent = TRUE,
        ensembl <- biomaRt::useMart("ensembl",
            dataset = "mmusculus_gene_ensembl",
            host    = "https://may2025.archive.ensembl.org")
    )
}
if (is.null(ensembl)) stop("Could not connect to Ensembl biomaRt after 5 attempts.")

all_ids <- biomaRt::getBM(attributes = "ensembl_gene_id", mart = ensembl)

genes <- NULL
attempts <- 0
while (is.null(genes) && attempts < 5) {
    attempts <- attempts + 1
    try(silent = TRUE,
        genes <- biomaRt::getBM(
            attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name",
                           "strand", "start_position", "end_position",
                           "gene_biotype", "go_id"),
            filters    = "ensembl_gene_id",
            values     = all_ids,
            mart       = ensembl,
            useCache   = TRUE)
    )
}
if (is.null(genes)) stop("biomaRt gene annotation query failed after 5 attempts.")

genes <- as_tibble(genes) %>%
    filter(!is.na(chromosome_name), chromosome_name %in% c(as.character(1:19), "X", "Y"))
cat(sprintf("  Retrieved %d rows for %d unique genes.\n",
    nrow(genes), n_distinct(genes$ensembl_gene_id)))

# ── 2. Gene coordinates file ──────────────────────────────────────────────────
cat("Writing gene coordinates...\n")
coords <- genes %>%
    filter(gene_biotype == "protein_coding") %>%
    distinct(chromosome_name, start_position, end_position,
             ensembl_gene_id, mgi_symbol) %>%
    drop_na()
write_tsv_nohdr(coords,
    file.path(args$outdir, "genes_coordinates_for_INRICH.txt"))
cat(sprintf("  %d protein-coding genes written.\n", nrow(coords)))

# ── 3. KEGG pathways ──────────────────────────────────────────────────────────
cat("Writing KEGG pathway gene sets...\n")
unique_ids <- unique(genes$ensembl_gene_id)
kegg_tbl <- AnnotationDbi::select(org.Mm.eg.db,
    columns  = "PATH",
    keys     = unique_ids,
    keytype  = "ENSEMBL") %>%
    as_tibble() %>%
    filter(!is.na(PATH)) %>%
    mutate(kdesc = paste0("KEGG map", PATH)) %>%
    dplyr::select(ENSEMBL, PATH, kdesc)
write_tsv_nohdr(kegg_tbl,
    file.path(args$outdir, "groups_KEGG_pathway_link_for_INRICH.txt"))
cat(sprintf("  %d gene-pathway pairs written.\n", nrow(kegg_tbl)))

# ── 4. GO terms (BP, CC, MF) ─────────────────────────────────────────────────
cat("Writing GO gene sets...\n")
go_tbl <- genes %>%
    distinct(ensembl_gene_id, go_id) %>%
    filter(go_id != "")

# Build GO-id → ontology + description lookup from GO.db
go_ids <- unique(go_tbl$go_id)
go_meta <- tibble(
    go_id = go_ids,
    ont   = vapply(go_ids, function(x) {
        g <- GO.db::GOTERM[[x]]
        if (!is.null(g)) g@Ontology else NA_character_
    }, character(1)),
    desc  = vapply(go_ids, function(x) {
        g <- GO.db::GOTERM[[x]]
        if (!is.null(g)) paste0(g@Term, " (", g@Ontology, ")") else NA_character_
    }, character(1))
)
go_tbl <- left_join(go_tbl, go_meta, by = "go_id") %>% filter(!is.na(ont))

for (ont in c("BP", "CC", "MF")) {
    out <- go_tbl %>%
        filter(ont == !!ont) %>%
        dplyr::select(ensembl_gene_id, go_id, desc)
    write_tsv_nohdr(out,
        file.path(args$outdir, paste0("groups_GO_", ont, "_terms_link_for_INRICH.txt")))
    cat(sprintf("  GO %s: %d gene-term pairs written.\n", ont, nrow(out)))
}

# ── 5. MGI → Ensembl mapping (shared by MP sections below) ───────────────────
# org.Mm.eg.db returns MGI IDs without the "MGI:" prefix (e.g. "1344407").
# Strip the prefix from extdata MGI columns at join time so both sides match.
mgi_map <- AnnotationDbi::select(org.Mm.eg.db,
    columns = "MGI",
    keys    = unique_ids,
    keytype = "ENSEMBL") %>%
    as_tibble() %>%
    filter(!is.na(MGI)) %>%
    mutate(MGI = sub(".*MGI:(\\d+).*", "\\1", MGI)) %>%
    filter(grepl("^\\d+$", MGI)) %>%
    distinct()

# ── 6. MP terms (Maayan lab gene-attribute edges) ─────────────────────────────
cat("Writing MP terms (Maayan lab)...\n")
mp_file <- file.path(args$extdata, "gene_attribute_edges.txt")
if (!file.exists(mp_file)) stop("Missing extdata file: ", mp_file)
# First row is a source/target header; second row is the actual column header.
mayphen <- read.table(mp_file, skip = 1, sep = "\t",
                      quote = "\"", header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(MGI.Accession = sub("^MGI:", "", MGI.Accession))
mp <- left_join(mayphen, mgi_map, by = c("MGI.Accession" = "MGI")) %>%
    dplyr::select(ENSEMBL, MPID, Phenotype) %>%
    filter(!is.na(ENSEMBL))
write_tsv_nohdr(mp,
    file.path(args$outdir, "groups_MP_terms_for_INRICH.txt"))
cat(sprintf("  %d gene-MP pairs written.\n", nrow(mp)))

# ── 7. MPMGI terms (MGI allele phenotype table) ───────────────────────────────
cat("Writing MPMGI terms...\n")
mpmgi_file <- file.path(args$extdata, "MPhenotype_MGenotype.csv")
if (!file.exists(mpmgi_file)) stop("Missing extdata file: ", mpmgi_file)
mphen <- read_csv(mpmgi_file,
    col_names = c("MGI", "name", "allele", "MPID", "MPname"),
    show_col_types = FALSE) %>%
    mutate(MGI = sub("^MGI:", "", MGI))
mp2 <- left_join(mphen, mgi_map, by = c("MGI" = "MGI")) %>%
    dplyr::select(ENSEMBL, MPID, MPname) %>%
    filter(!is.na(ENSEMBL))
write_tsv_nohdr(mp2,
    file.path(args$outdir, "groups_MPMGI_terms_for_INRICH.txt"))
cat(sprintf("  %d gene-MPMGI pairs written.\n", nrow(mp2)))

# ── 8. Brain expression groups ────────────────────────────────────────────────
cat("Writing brain expression groups...\n")
brain_file <- file.path(args$extdata,
    "groups_from_41598_2016_BFsrep19274_MOESM4_ESM.txt")
if (!file.exists(brain_file)) stop("Missing extdata file: ", brain_file)
brain_wide <- read.delim(brain_file, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE)
brain_long <- brain_wide %>%
    pivot_longer(-V1, values_to = "mgi_symbol") %>%
    dplyr::select(Group = V1, mgi_symbol) %>%
    filter(mgi_symbol != "") %>%
    mutate(grname = sapply(Group, function(x)
        paste(strsplit(x, " ")[[1]][1:2], collapse = "_")))

# Map mgi_symbol → ensembl_gene_id via the biomaRt result
sym_map <- genes %>% distinct(ensembl_gene_id, mgi_symbol) %>% filter(mgi_symbol != "")
brain_out <- left_join(brain_long, sym_map, by = "mgi_symbol",
                       relationship = "many-to-many") %>%
    filter(!is.na(ensembl_gene_id)) %>%
    dplyr::select(ensembl_gene_id, grname, Group) %>%
    distinct()
write_tsv_nohdr(brain_out,
    file.path(args$outdir, "groups_brain_expression_link_for_INRICH.txt"))
cat(sprintf("  %d gene-group pairs written across %d groups.\n",
    nrow(brain_out), n_distinct(brain_out$Group)))

# ── 9. SFARI ASD gene sets ────────────────────────────────────────────────────
if (!args$no_sfari) {
    cat("Writing SFARI gene sets (requires human biomaRt for ortholog mapping)...\n")
    sfari_file <- file.path(args$extdata,
        "SFARI-Gene_genes_01-13-2021release_01-22-2021export.csv")
    if (!file.exists(sfari_file)) {
        cat("  SFARI file not found, skipping.\n")
    } else {
        sft <- read_csv(sfari_file, show_col_types = FALSE) %>%
            filter(!is.na(`ensembl-id`))

        human_mart <- NULL
        mouse_mart <- NULL
        try(silent = TRUE, {
            human_mart <- biomaRt::useMart("ensembl",
                dataset = "hsapiens_gene_ensembl",
                host    = "https://may2025.archive.ensembl.org")
            mouse_mart <- biomaRt::useMart("ensembl",
                dataset = "mmusculus_gene_ensembl",
                host    = "https://may2025.archive.ensembl.org")
        })
        if (is.null(human_mart) || is.null(mouse_mart)) {
            cat("  Could not connect to human biomaRt, skipping SFARI.\n")
        } else {
            orthologs <- tryCatch(
                biomaRt::getLDS(
                    attributes  = "ensembl_gene_id",
                    filters     = "ensembl_gene_id",
                    values      = sft$`ensembl-id`,
                    mart        = human_mart,
                    attributesL = "ensembl_gene_id",
                    martL       = mouse_mart,
                    uniqueRows  = TRUE),
                error = function(e) {
                    cat(sprintf("  getLDS failed (%s), skipping SFARI.\n",
                                conditionMessage(e)))
                    NULL
                })
            if (is.null(orthologs)) {
                cat("  SFARI skipped due to getLDS error.\n")
            } else {
            colnames(orthologs) <- c("human_ensembl", "mouse_ensembl")
            sft2 <- left_join(sft, orthologs,
                              by = c("ensembl-id" = "human_ensembl")) %>%
                filter(!is.na(mouse_ensembl))
            sfari_sets <- bind_rows(
                tibble(gene  = sft2$mouse_ensembl,
                       trait = "SFARI_all",
                       desc  = "SFARI ASD 1-3 gene list"),
                tibble(gene  = sft2$mouse_ensembl[
                                    !is.na(sft2$`gene-score`) &
                                    sft2$`gene-score` <= 2],
                       trait = "SFARI_1_2",
                       desc  = "SFARI ASD 1-2 gene list"),
                tibble(gene  = sft2$mouse_ensembl[
                                    !is.na(sft2$`gene-score`) &
                                    sft2$`gene-score` == 1],
                       trait = "SFARI_1",
                       desc  = "SFARI ASD 1 gene list")
            )
            write_tsv_nohdr(sfari_sets,
                file.path(args$outdir, "groups_SFARI_terms_for_INRICH.txt"))
            cat(sprintf("  %d gene-set rows written (%d unique mouse genes).\n",
                nrow(sfari_sets), n_distinct(sfari_sets$gene)))
            } # end if (!is.null(orthologs))
        }
    }
} else {
    cat("Skipping SFARI (--no-sfari).\n")
}

cat(sprintf("\nDone. Static INRICH files written to: %s\n", args$outdir))
