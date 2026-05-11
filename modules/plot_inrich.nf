process PLOT_INRICH {
    tag "${name}"
    label 'low'

    publishDir "${params.outdir}/inrich/plots", mode: 'copy'

    container params.container

    input:
    val  name
    path result_files   // all *.out.inrich collected from RUN_INRICH
    path yaml

    output:
    path "*.pdf", emit: plots, optional: true

    script:
    """
    Rscript ${projectDir}/src/plot_inrich.R \\
        --files  ${result_files.join(' ')} \\
        --yaml   ${yaml}                   \\
        --name   ${name}                   \\
        --outdir .                         \\
        --minpval ${params.inrich_minpval}
    """
}
