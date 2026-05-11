#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// Boliwood GWAS pipeline
//
// Stages (added one by one):
//   1. GEMMA kinship  — compute LOCO kinship matrices
//   2. GEMMA LMM      — run LMM per trait per chromosome, concat results
//   3. SuSiE          — fine-map each significant locus per trait
//   4. Manhattan plots — per-trait plots with credible sets overlaid
//   5. Group Manhattan — ACAT-combined p-values per phenotype group
//   6. Pleiotropy matrix — cross-trait p-value matrix at SuSiE loci
//   7. Peak clustering  — cluster loci by cross-trait p-value patterns
// ─────────────────────────────────────────────────────────────────────────────

// ── Module imports (uncommented as each stage is added) ─────────────────────
include { PREPARE_INPUT         } from './modules/prepare_input'
include { KINSHIP_FULL          } from './modules/kinship_full'
include { KINSHIP_LOCO          } from './modules/kinship_loco'
include { HERITABILITY              } from './modules/heritability'
include { COLLECT_HERITABILITY      } from './modules/heritability'
include { KINSHIP_LDAK              } from './modules/kinship_ldak'
include { HERITABILITY_LDAK         } from './modules/heritability_ldak'
include { KINSHIP_GCTA              } from './modules/kinship_gcta'
include { HERITABILITY_GCTA         } from './modules/heritability_gcta'
include { COLLECT_HERITABILITY_GCTA } from './modules/heritability_gcta'
include { GEMMA_LMM             } from './modules/gemma_lmm'
include { COLLECT_GEMMA_LMM     } from './modules/gemma_lmm'
include { PERMUTATION           } from './modules/permutation'
include { CALC_PERM_THRESHOLD   } from './modules/permutation'
include { CLUMP                 } from './modules/clump'
include { CLUMP_COMBINED        } from './modules/clump_combined'
include { SUSIE_FINEMAP         } from './modules/susie_finemap'
include { GET_GENES             } from './modules/get_genes'
include { COMBINE_GWAS          } from './modules/combine_gwas'
include { PLOT_MANHATTAN        } from './modules/plot_manhattan'
include { GENETIC_CORR          } from './modules/genetic_corr'
include { COLLECT_GENETIC_CORR  } from './modules/genetic_corr'
// include { PLOT_GROUP_MANHATTAN } from './modules/plot_group_manhattan'
// include { PLEIOTROPY_PVALMAT   } from './modules/pleiotropy_pvalmat'
include { CLUSTER_PEAKS         } from './modules/cluster_peaks'
include { PLOT_PVE              } from './modules/plot_pve'
include { PLOT_LD_DECAY         } from './modules/plot_ld_decay'
include { WRITE_INRICH_SNP_MAP
          WRITE_INRICH_INTERVALS
          WRITE_INRICH_INTERVALS_COMBINED
          RUN_INRICH              } from './modules/inrich'
include { PLOT_INRICH           } from './modules/plot_inrich'

// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // ── Input validation ────────────────────────────────────────────────────
    if ( !params.input      ) error "Please provide --input (phenotype CSV)"
    if ( !params.yaml       ) error "Please provide --yaml (run config YAML)"
    if ( !params.mda && !params.genotypes ) error "Please provide --genotypes (HMDP genotype CSV) or set --mda true"
    if ( !params.name       ) error "Please provide --name (run name)"

    // ── Stage 0: Prepare inputs ─────────────────────────────────────────────
    def geno_ch = params.mda
        ? Channel.fromPath("${params.extdata}/snp_*.csv.gz").collect()
        : Channel.value(file(params.genotypes))

    PREPARE_INPUT(
        file(params.input),
        file(params.yaml),
        geno_ch,
        params.name
    )

    // Derive trait indices from actual phenotypes_order.txt output rather than
    // params.n_traits, so dropped traits (all-NA/zero after scaling) don't
    // cause off-by-one errors in downstream array jobs.
    trait_idx_ch = PREPARE_INPUT.out.pheno_order
        .flatMap { f ->
            def n = f.readLines().drop(1)*.trim().findAll { it }.size()
            (1..n).toList()
        }

    // ── Stage 0b: LD decay plot ─────────────────────────────────────────────
    PLOT_LD_DECAY(
        PREPARE_INPUT.out.plink,
        file("${projectDir}/src/plot_ld_decay.R"),
        file("${projectDir}/src/aesthetics.R"),
        file("${projectDir}/src/postprocess_setup.R")
    )

    // ── Stage 1: GEMMA kinship ──────────────────────────────────────────────
    KINSHIP_FULL(
        PREPARE_INPUT.out.plink,
        PREPARE_INPUT.out.pheno
    )

    KINSHIP_LOCO(
        PREPARE_INPUT.out.plink.combine(Channel.of(1..19, 'X')),
        PREPARE_INPUT.out.pheno
    )

    // ── Stage 2: Heritability (GEMMA) ───────────────────────────────────────
    HERITABILITY(
        PREPARE_INPUT.out.plink.combine(trait_idx_ch),
        PREPARE_INPUT.out.pheno,
        PREPARE_INPUT.out.covars,
        KINSHIP_FULL.out.kinship
    )

    COLLECT_HERITABILITY(
        params.name,
        HERITABILITY.out.tsv.collect(),
        PREPARE_INPUT.out.pheno_order
    )

    // ── Stage 2: Heritability (LDAK) ────────────────────────────────────────
    KINSHIP_LDAK(
        PREPARE_INPUT.out.plink
    )

    HERITABILITY_LDAK(
        params.name,
        PREPARE_INPUT.out.pheno_plink,
        PREPARE_INPUT.out.covars_plink,
        KINSHIP_LDAK.out.kinship,
        PREPARE_INPUT.out.pheno_order
    )

    // ── Stage 2b: Heritability (GCTA) ───────────────────────────────────────
    KINSHIP_GCTA(
        PREPARE_INPUT.out.plink
    )

    HERITABILITY_GCTA(
        trait_idx_ch,
        KINSHIP_GCTA.out.kinship,
        PREPARE_INPUT.out.pheno_plink,
        PREPARE_INPUT.out.covars_plink
    )

    COLLECT_HERITABILITY_GCTA(
        params.name,
        HERITABILITY_GCTA.out.csv.collect(),
        PREPARE_INPUT.out.pheno_order
    )

    // ── Stage 2c: Heritability plots ────────────────────────────────────────
    PLOT_PVE(
        COLLECT_HERITABILITY.out.heritability,
        COLLECT_HERITABILITY_GCTA.out.heritability,
        HERITABILITY_LDAK.out.heritability,
        params.name,
        file(params.yaml),
        file("${projectDir}/src/plot_pve.R"),
        file("${projectDir}/src/aesthetics.R"),
        file("${projectDir}/src/postprocess_setup.R")
    )

    // ── Stage 3: Permutation threshold ─────────────────────────────────────
    PERMUTATION(
        PREPARE_INPUT.out.plink.combine(Channel.of(1..params.n_perms)),
        PREPARE_INPUT.out.pheno,
        PREPARE_INPUT.out.covars,
        KINSHIP_FULL.out.kinship,
        params.perm_trait,
        params.pval_type
    )

    CALC_PERM_THRESHOLD(
        PERMUTATION.out.minp.collect()
    )

    // ── Stage 4: GEMMA LMM ─────────────────────────────────────────────────
    // Join chr_bfiles and kinship on chr, then combine with trait indices
    lmm_input_ch = trait_idx_ch
        .combine(
            KINSHIP_LOCO.out.chr_bfiles.join(KINSHIP_LOCO.out.kinship)
        )
        .map { trait_idx, chr, bed, bim, fam, kinship ->
            tuple(trait_idx, chr, bed, bim, fam, kinship)
        }

    GEMMA_LMM(
        lmm_input_ch,
        PREPARE_INPUT.out.pheno,
        PREPARE_INPUT.out.covars
    )

    // Group per-chr results by trait, then concatenate
    COLLECT_GEMMA_LMM(
        GEMMA_LMM.out.assoc.groupTuple()
    )

    // ── Stage 4b: Combine GWAS results into wide-format matrices ────────────
    COMBINE_GWAS(
        params.name,
        COLLECT_GEMMA_LMM.out.assoc.map { trait_idx, f -> f }.collect(),
        PREPARE_INPUT.out.pheno_order
    )

    // ── Stage 5: LD clumping ────────────────────────────────────────────────
    CLUMP(
        COLLECT_GEMMA_LMM.out.assoc,
        PREPARE_INPUT.out.plink,
        CALC_PERM_THRESHOLD.out.threshold,
        params.pval_type
    )

    // ── Stage 5b: Combined LD clumping (all-traits + per-group) ─────────────
    // Identifies pleiotropic QTLs and per-group QTL counts.
    // Runs once on the consolidated p-value matrix; no extra R script needed.
    CLUMP_COMBINED(
        COMBINE_GWAS.out.pvals,
        PREPARE_INPUT.out.trait_groups,
        PREPARE_INPUT.out.plink,
        CALC_PERM_THRESHOLD.out.threshold,
        params.pval_type
    )

    // ── Stage 6: SuSiE fine-mapping ─────────────────────────────────────────
    // When run_susie = false, skip SUSIE_FINEMAP and route CLUMP outputs downstream.
    def inrich_mode = params.run_susie ? "susie" : "clump"

    if (params.run_susie) {
        susie_input_ch = CLUMP.out.clumped
            .join(COLLECT_GEMMA_LMM.out.assoc)

        SUSIE_FINEMAP(
            susie_input_ch,
            PREPARE_INPUT.out.plink,
            PREPARE_INPUT.out.pheno,
            CALC_PERM_THRESHOLD.out.threshold
        )

        // Flat list of all credible-set files for CLUSTER_PEAKS / WRITE_INRICH_INTERVALS_COMBINED
        susie_files_ch  = SUSIE_FINEMAP.out.credible_sets
                              .map { trait_idx, files -> files }
                              .flatten()
                              .collect()
        // Per-trait channel for WRITE_INRICH_INTERVALS
        inrich_locus_ch = SUSIE_FINEMAP.out.credible_sets
    } else {
        // Use per-trait clump files as the locus source for all downstream steps.
        // Staged in susie_inputs/ inside CLUSTER_PEAKS to avoid filename collisions.
        susie_files_ch  = CLUMP.out.clumped.map { trait_idx, f -> f }.collect()
        // Join clump and assoc on trait_idx so both files are staged together in
        // WRITE_INRICH_INTERVALS — the R script discovers the assoc file by glob.
        inrich_locus_ch = CLUMP.out.clumped
            .join(COLLECT_GEMMA_LMM.out.assoc)
            .map { trait_idx, clumped, assoc -> tuple(trait_idx, [clumped, assoc]) }
    }

    // ── Stage 7: Gene annotation ─────────────────────────────────────────────
    // When run_susie = false, always fall back to clump files regardless of annot_input param.
    def annot_mode
    if (params.run_susie && params.annot_input == "susie") {
        annot_mode     = "susie"
        annot_input_ch = SUSIE_FINEMAP.out.credible_sets.map { trait_idx, files -> files }.collect()
    } else {
        annot_mode     = "clump"
        annot_input_ch = CLUMP.out.clumped.map { trait_idx, file -> file }.collect()
    }

    GET_GENES(annot_input_ch, annot_mode)

    // ── Stage 8: Genetic correlations (LDAK bivariate REML) ─────────────────
    // Build lower-triangular pairs from corr_traits.txt subset (or all traits)
    def corr_trait_names = params.corr_traits
        ? file(params.corr_traits).readLines().findAll { it.trim() && !it.startsWith('#') }
        : null

    corr_pairs_ch = PREPARE_INPUT.out.pheno_order
        .map { f ->
            def all_traits = f.readLines().drop(1)*.trim().findAll { it }
            def subset     = corr_trait_names ?: all_traits
            def name_to_idx = [:]
            all_traits.eachWithIndex { t, idx -> name_to_idx[t] = idx + 1 }
            def indices = subset.collect { name_to_idx[it] }.findAll { it != null }
            def missing = subset.findAll { !name_to_idx.containsKey(it) }
            if (missing) log.warn "corr_traits: ${missing.size()} name(s) not found in phenotypes_order.txt: ${missing}"
            // Lower-triangular: i > j
            def pairs = []
            indices.eachWithIndex { a, ai ->
                indices[0..<ai].each { b -> pairs << [a, b] }
            }
            pairs
        }
        .flatMap { it }
        .map { pair -> tuple(pair[0], pair[1]) }

    GENETIC_CORR(
        corr_pairs_ch,
        PREPARE_INPUT.out.pheno_plink,
        PREPARE_INPUT.out.covars_plink,
        KINSHIP_GCTA.out.kinship
    )

    COLLECT_GENETIC_CORR(
        params.name,
        GENETIC_CORR.out.hsq.collect(),
        PREPARE_INPUT.out.pheno_order
    )

    // ── Stage 9: Manhattan plots (individual + group) ────────────────────────
    PLOT_MANHATTAN(
        params.name,
        COMBINE_GWAS.out.pvals,
        file(params.yaml),
        CALC_PERM_THRESHOLD.out.threshold,
        params.pval_thresh,
        params.pval_type,
        file("${projectDir}/src/plot_manhattan.R"),
        file("${projectDir}/src/aesthetics.R"),
        file("${projectDir}/src/postprocess_setup.R"),
        GET_GENES.out.annotations
    )

    // ── Stage 11: INRICH enrichment analysis ─────────────────────────────────
    WRITE_INRICH_SNP_MAP(
        PREPARE_INPUT.out.plink
    )

    // Per-trait interval files — locus source is SuSiE credible sets or CLUMP depending on run_susie
    WRITE_INRICH_INTERVALS(
        inrich_locus_ch,
        PREPARE_INPUT.out.pheno_order,
        inrich_mode,
        PREPARE_INPUT.out.plink
    )

    // Stage all pre-generated gene-set files from inrich_static/
    inrich_static_ch = Channel.fromPath("${params.inrich_static}/*.txt").collect()

    // Per-trait intervals → (label, file) tuples
    trait_inrich_ch = WRITE_INRICH_INTERVALS.out.intervals
        .flatten()
        .map { f ->
            def label = f.name.replaceAll('^intervals_', '').replaceAll('_for_INRICH\\.txt$', '')
            tuple(label, f)
        }

    // ── Stage 10: Peak clustering ─────────────────────────────────────────────
    CLUSTER_PEAKS(
        params.name,
        COMBINE_GWAS.out.pvals,
        susie_files_ch,
        file(params.yaml),
        file("${projectDir}/src/cluster_peaks.R"),
        file("${projectDir}/src/aesthetics.R"),
        file("${projectDir}/src/postprocess_setup.R"),
        CLUMP_COMBINED.out.clumped,
        CLUMP.out.clumped.map { trait_idx, f -> f }.collect(),
        params.pval_type,
        params.run_susie
    )

    // ── Stage 11 (cont.): per-group + per-cluster INRICH (after clustering) ──
    WRITE_INRICH_INTERVALS_COMBINED(
        susie_files_ch,
        PREPARE_INPUT.out.pheno_order,
        PREPARE_INPUT.out.trait_groups,
        CLUSTER_PEAKS.out.csv.collect(),
        inrich_mode,
        COLLECT_GEMMA_LMM.out.assoc.map { trait_idx, f -> f }.collect(),
        PREPARE_INPUT.out.plink
    )

    combined_inrich_ch = WRITE_INRICH_INTERVALS_COMBINED.out.intervals
        .flatten()
        .map { f ->
            def label = f.name.replaceAll('^intervals_', '').replaceAll('_for_INRICH\\.txt$', '')
            tuple(label, f)
        }

    RUN_INRICH(
        trait_inrich_ch.mix(combined_inrich_ch),
        WRITE_INRICH_SNP_MAP.out.snp_map,
        inrich_static_ch
    )

    // ── Stage 12: Plot INRICH results ─────────────────────────────────────────
    PLOT_INRICH(
        params.name,
        RUN_INRICH.out.results.collect(),
        file(params.yaml)
    )
}
