process PREPARE_INPUT {
    tag "${name}"
    label 'medium'

    publishDir "${params.outdir}/input", mode: 'copy'

    container params.r_container

    input:
    path  csv
    path  yaml
    path  genotypes
    val   name

    output:
    tuple val(name), path("${name}.bed"), path("${name}.bim"), path("${name}.fam"),  emit: plink
    path "pheno_${name}.txt",     emit: pheno
    path "anno_${name}.txt",      emit: anno
    path "covars_${name}.txt",    emit: covars
    path "annotations.csv",       emit: annotations
    path "phenotypes_order.txt",  emit: pheno_order
    path "raw_phenotypes.csv",    emit: raw_pheno

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
        --missing    ${params.missing}
    """
}
