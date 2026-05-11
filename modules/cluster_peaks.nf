process CLUSTER_PEAKS {
    tag "${name}"
    label 'low'

    publishDir "${params.outdir}/plots/clusters", mode: 'copy'

    container params.container

    input:
    val  name
    path pvals
    path susie_files, stageAs: 'susie_inputs/*'   // staged in sub-dir to avoid clump filename collisions
    path yaml
    path script        // cluster_peaks.R
    path aesthetics    // aesthetics.R — staged alongside script so dirname() lookup works
    path setup         // postprocess_setup.R — staged alongside script
    path clumped_files // *.clumped from CLUMP_COMBINED (all.clumped + per-group)
    path clump_files   // trait<N>.clumped from CLUMP (per-trait, collected)
    val  pval_type     // which p-value column for clustering: p_lrt, p_score, or p_wald
    val  run_susie     // false = use PLINK clump leads instead of SuSiE credible sets

    output:
    path "*.pdf",  emit: plots
    path "*.svg",  emit: svg,     optional: true
    path "*.csv",  emit: csv

    script:
    def susie_arg = run_susie
        ? "--susie ${susie_files.join(' ')}"
        : "--no_susie"
    """
    Rscript ${script} \\
        --pvals       ${pvals}                     \\
        ${susie_arg}                               \\
        --name        ${name}                      \\
        --yaml        ${yaml}                      \\
        --outdir      .                            \\
        --clusters    ${params.k_clusters}         \\
        --pvalthr     ${params.pvalthr_clusters}   \\
        --clumped_dir .                            \\
        --clump_dir   .                            \\
        --pval_type   ${pval_type}
    """
}
