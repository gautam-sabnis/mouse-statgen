process GET_GENES {
    tag "gene_annotation"
    label 'medium'

    publishDir "${params.outdir}/annotations", mode: 'copy'

    container params.container

    input:
    path input_files   // either all SuSiE locus files or all clump files

    output:
    path "gene_annotations.csv", emit: annotations

    script:
    """
    Rscript ${projectDir}/src/get_genes.R \\
        --mode    ${params.annot_input} \\
        --input   ${input_files} \\
        --dist    ${params.locus_window} \\
        --outdir  .
    """
}
