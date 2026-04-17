process CLUSTER_PEAKS {
    tag "${name}"
    label 'low'

    publishDir "${params.outdir}/plots/clusters", mode: 'copy'

    container params.container

    input:
    val  name
    path pvals
    path susie_files   // trait*_locus_*.txt from SUSIE_FINEMAP (collected)
    path yaml
    path script        // cluster_peaks.R
    path aesthetics    // aesthetics.R — staged alongside script so dirname() lookup works
    path setup         // postprocess_setup.R — staged alongside script
    path clumped_files // *.clumped from CLUMP_COMBINED (all.clumped + per-group)
    path clump_files   // trait<N>.clumped from CLUMP (per-trait, collected)

    output:
    path "*.pdf",  emit: plots
    path "*.svg",  emit: svg,     optional: true
    path "*.csv",  emit: csv

    script:
    """
    Rscript ${script} \\
        --pvals       ${pvals}                     \\
        --susie       ${susie_files.join(' ')}     \\
        --name        ${name}                      \\
        --yaml        ${yaml}                      \\
        --outdir      .                            \\
        --clusters    ${params.k_clusters}         \\
        --pvalthr     ${params.pvalthr_clusters}   \\
        --clumped_dir .                            \\
        --clump_dir   .
    """
}
