process KINSHIP_LDAK {
    tag "${name}"
    label 'medium'

    publishDir "${params.outdir}/kinship", mode: 'copy'

    container params.container

    input:
    tuple val(name), path(bed), path(bim), path(fam)

    output:
    tuple path("ldak_kinship.grm.bin"),
          path("ldak_kinship.grm.id"),
          path("ldak_kinship.grm.details"), emit: kinship

    script:
    """
    ldak --calc-kins-direct ldak_kinship \\
         --bfile ${name} \\
         --power -0.25
    """
}
