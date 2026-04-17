process SUSIE_FINEMAP {
    tag "trait${trait_idx}"
    label 'medium'

    publishDir "${params.outdir}/susie", mode: 'copy'

    container params.container

    input:
    tuple val(trait_idx), path(clumped), path(assoc)
    tuple val(name), path(bed), path(bim), path(fam)
    path pheno
    path threshold

    output:
    tuple val(trait_idx), path("trait${trait_idx}_locus_*.txt"), optional: true, emit: credible_sets

    script:
    """
    Rscript ${projectDir}/src/susie_finemap.R \\
        --bfile     ${name} \\
        --pheno     ${pheno} \\
        --trait     ${trait_idx} \\
        --gwas      ${assoc} \\
        --clump     ${clumped} \\
        --threshold ${threshold} \\
        --window    ${params.locus_window} \\
        --L         ${params.susie_L} \\
        --outdir    .
    """
}
