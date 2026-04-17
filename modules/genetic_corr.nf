process GENETIC_CORR {
    tag "rg_${i}_${j}"
    label 'low'

    publishDir "${params.outdir}/genetic_corr/pairs", mode: 'copy'

    container params.container

    input:
    tuple val(i), val(j)
    path  pheno_plink
    path  covars_plink
    tuple path(grm_bin), path(grm_n_bin), path(grm_id)

    output:
    path "rg_${i}_${j}.hsq", emit: hsq, optional: true

    script:
    """
    gcta --reml-bivar ${i} ${j} \\
         --grm gcta_kinship \\
         --pheno ${pheno_plink} \\
         --qcovar ${covars_plink} \\
         --reml-bivar-lrt-rg 0 \\
         --out rg_${i}_${j} \\
         --thread-num ${task.cpus}
    """
}

process COLLECT_GENETIC_CORR {
    tag "${name}"
    label 'low'

    publishDir "${params.outdir}/genetic_corr", mode: 'copy'

    container params.container

    input:
    val  name
    path hsq_files, stageAs: "hsq/*"
    path pheno_order

    output:
    path "genetic_corr_${name}.csv", emit: corr_matrix
    path "genetic_pval_${name}.csv", emit: pval_matrix

    script:
    """
    Rscript ${projectDir}/src/collect_genetic_corr.R \\
        --hsq_dir     hsq \\
        --pheno_order ${pheno_order} \\
        --name        ${name} \\
        --outdir      .
    """
}
