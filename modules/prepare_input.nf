process PREPARE_INPUT {
    tag "${name}"
    label 'high'

    publishDir "${params.outdir}/input", mode: 'copy'

    container params.container

    input:
    path  csv
    path  yaml
    path  genotypes
    val   name

    output:
    tuple val(name), path("${name}.bed"), path("${name}.bim"), path("${name}.fam"),  emit: plink
    path "pheno_${name}.txt",          emit: pheno
    path "pheno_plink_${name}.txt",    emit: pheno_plink
    path "anno_${name}.txt",           emit: anno
    path "covars_${name}.txt",         emit: covars
    path "covars_plink_${name}.txt",   emit: covars_plink
    path "annotations.csv",            emit: annotations
    path "phenotypes_order.txt",       emit: pheno_order
    path "raw_phenotypes.csv",         emit: raw_pheno
    path "trait_groups_${name}.tsv",   emit: trait_groups

    script:
    """
    Rscript ${projectDir}/src/prepare_and_make_plink.R \\
        --input      ${csv}        \\
        --yaml       ${yaml}       \\
        --genotypes  ${genotypes}  \\
        --name       ${name}       \\
        --outdir     .             \\
        --downsample ${params.downsample} \\
        --MAF        ${params.maf}        \\
        --missing    ${params.missing}    \\
        ${params.qqnorm ? '--qqnorm' : ''}
    """
}
