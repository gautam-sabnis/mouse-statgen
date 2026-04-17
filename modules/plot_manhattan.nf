process PLOT_MANHATTAN {
    tag "${name}"
    label 'low'

    publishDir "${params.outdir}/plots/manhattan", mode: 'copy'

    container params.container

    input:
    val  name
    path pvals
    path yaml
    path threshold   // single-value file from CALC_PERM_THRESHOLD; pass 'NO_FILE' to use fallback
    val  pval_thresh // fallback threshold from params (used when threshold file is absent)
    path script      // plot_manhattan.R — declared so cache is invalidated on script changes
    path aesthetics  // aesthetics.R — staged alongside script so dirname() lookup works
    path setup       // postprocess_setup.R — staged alongside script
    path annot_file  // gene_annotations.csv from GET_GENES; pass 'NO_FILE' to skip gene labels

    output:
    path "*.pdf", emit: plots

    script:
    def thresh_arg = (threshold.name != 'NO_FILE')
        ? "--threshold_file ${threshold}"
        : "--threshold ${pval_thresh}"
    def annot_arg = (annot_file.name != 'NO_FILE')
        ? "--annot_file ${annot_file} --label_gap ${params.label_gap}"
        : ""
    """
    Rscript ${script} \\
        --pvals   ${pvals}  \\
        --yaml    ${yaml}   \\
        --name    ${name}   \\
        --outdir  .         \\
        ${thresh_arg}       \\
        ${annot_arg}
    """
}
