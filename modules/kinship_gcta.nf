process KINSHIP_GCTA {
    tag "${name}"
    label 'medium'

    publishDir "${params.outdir}/kinship", mode: 'copy'

    container params.container

    input:
    tuple val(name), path(bed), path(bim), path(fam)

    output:
    tuple path("gcta_kinship.grm.bin"),
          path("gcta_kinship.grm.N.bin"),
          path("gcta_kinship.grm.id"), emit: kinship

    script:
    """
    gcta --bfile ${name} \\
         --make-grm \\
         --out gcta_kinship \\
         --thread-num ${task.cpus}
    """
}
