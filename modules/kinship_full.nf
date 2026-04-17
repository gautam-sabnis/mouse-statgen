process KINSHIP_FULL {
    tag "${name}"
    label 'medium'

    publishDir "${params.outdir}/kinship", mode: 'copy'

    container params.container

    input:
    tuple val(name), path(bed), path(bim), path(fam)
    path pheno

    output:
    path "kin_full.cXX.txt", emit: kinship

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}

    gemma -bfile ${name} \\
          -p ${pheno} \\
          -n 1 \\
          -gk 1 \\
          -o kin_full
    mv output/*.cXX.txt .
    """
}
