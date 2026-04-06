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
include { PREPARE_INPUT  } from './modules/prepare_input'
// include { GEMMA_KINSHIP } from './modules/gemma_kinship'
// include { GEMMA_LMM     } from './modules/gemma_lmm'
// include { SUSIE_FINEMAP } from './modules/susie_finemap'
// include { PLOT_MANHATTAN       } from './modules/plot_manhattan'
// include { PLOT_GROUP_MANHATTAN } from './modules/plot_group_manhattan'
// include { PLEIOTROPY_PVALMAT   } from './modules/pleiotropy_pvalmat'
// include { CLUSTER_PEAKS        } from './modules/cluster_peaks'

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
    // (add here — takes PREPARE_INPUT.out.plink)

    // ── Stage 2: GEMMA LMM ─────────────────────────────────────────────────
    // (add here)

    // ── Stage 3: SuSiE fine-mapping ─────────────────────────────────────────
    // (add here)

    // ── Stage 4: Per-trait Manhattan plots ──────────────────────────────────
    // (add here)

    // ── Stage 5: Group Manhattan plots (ACAT) ───────────────────────────────
    // (add here)

    // ── Stage 6: Cross-trait pleiotropy matrix ───────────────────────────────
    // (add here)

    // ── Stage 7: Peak clustering ─────────────────────────────────────────────
    // (add here)
}
