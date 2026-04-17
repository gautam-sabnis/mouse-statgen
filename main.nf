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

// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // ── Input validation ────────────────────────────────────────────────────
    if ( !params.input      ) error "Please provide --input (phenotype CSV)"
    if ( !params.yaml       ) error "Please provide --yaml (run config YAML)"
    if ( !params.genotypes  ) error "Please provide --genotypes (HMDP genotype CSV)"
    if ( !params.name       ) error "Please provide --name (run name)"

    // ── Stage 0: Prepare inputs ─────────────────────────────────────────────
    PREPARE_INPUT(
        file(params.input),
        file(params.yaml),
        file(params.genotypes),
        params.name
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
        PREPARE_INPUT.out.plink.combine(Channel.of(1..params.n_traits)),
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
        Channel.of(1..params.n_traits),
        KINSHIP_GCTA.out.kinship,
        PREPARE_INPUT.out.pheno_plink,
        PREPARE_INPUT.out.covars_plink
    )

    COLLECT_HERITABILITY_GCTA(
        params.name,
        HERITABILITY_GCTA.out.csv.collect(),
        PREPARE_INPUT.out.pheno_order
    )

    // ── Stage 3: Permutation threshold ─────────────────────────────────────
    PERMUTATION(
        PREPARE_INPUT.out.plink.combine(Channel.of(1..params.n_perms)),
        PREPARE_INPUT.out.pheno,
        PREPARE_INPUT.out.covars,
        KINSHIP_FULL.out.kinship,
        params.perm_trait
    )

    CALC_PERM_THRESHOLD(
        PERMUTATION.out.minp.collect()
    )

    // ── Stage 4: GEMMA LMM ─────────────────────────────────────────────────
    // Join chr_bfiles and kinship on chr, then combine with trait indices
    // Produces 69 × 19 = 1311 parallel jobs
    lmm_input_ch = Channel.of(1..params.n_traits)
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
        CALC_PERM_THRESHOLD.out.threshold
    )

    // ── Stage 5b: Combined LD clumping (all-traits + per-group) ─────────────
    // Identifies pleiotropic QTLs and per-group QTL counts.
    // Runs once on the consolidated p-value matrix; no extra R script needed.
    CLUMP_COMBINED(
        COMBINE_GWAS.out.pvals,
        PREPARE_INPUT.out.trait_groups,
        PREPARE_INPUT.out.plink,
        CALC_PERM_THRESHOLD.out.threshold
    )

    // ── Stage 6: SuSiE fine-mapping ─────────────────────────────────────────
    // Join clump and assoc channels on trait_idx, use permutation threshold
    susie_input_ch = CLUMP.out.clumped
        .join(COLLECT_GEMMA_LMM.out.assoc)

    SUSIE_FINEMAP(
        susie_input_ch,
        PREPARE_INPUT.out.plink,
        PREPARE_INPUT.out.pheno,
        CALC_PERM_THRESHOLD.out.threshold
    )

    // ── Stage 7: Gene annotation ─────────────────────────────────────────────
    if (params.annot_input == "susie") {
        annot_input_ch = SUSIE_FINEMAP.out.credible_sets.map { trait_idx, files -> files }.collect()
    } else {
        annot_input_ch = CLUMP.out.clumped.map { trait_idx, file -> file }.collect()
    }

    GET_GENES(annot_input_ch)

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
        file("${projectDir}/src/plot_manhattan.R"),
        file("${projectDir}/src/aesthetics.R"),
        file("${projectDir}/src/postprocess_setup.R"),
        GET_GENES.out.annotations
    )

    // ── Stage 10: Peak clustering ─────────────────────────────────────────────
    CLUSTER_PEAKS(
        params.name,
        COMBINE_GWAS.out.pvals,
        SUSIE_FINEMAP.out.credible_sets
            .map { trait_idx, files -> files }
            .flatten()
            .collect(),
        file(params.yaml),
        file("${projectDir}/src/cluster_peaks.R"),
        file("${projectDir}/src/aesthetics.R"),
        file("${projectDir}/src/postprocess_setup.R"),
        CLUMP_COMBINED.out.clumped,
        CLUMP.out.clumped.map { trait_idx, f -> f }.collect()
    )
}
