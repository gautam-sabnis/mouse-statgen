process GET_GENES {
    tag "gene_annotation"
    label 'medium'

    publishDir "${params.outdir}/annotations", mode: 'copy'

    container params.container

    input:
    path input_files   // either all SuSiE locus files or all clump files
    val  mode          // "susie" or "clump" — resolved by main.nf, not taken from params directly

    output:
    path "gene_annotations.csv", emit: annotations

    script:
    """
    Rscript ${projectDir}/src/get_genes.R \\
        --mode    ${mode} \\
        --input   ${input_files} \\
        --dist    ${params.locus_window} \\
        --outdir  .
    """
}
