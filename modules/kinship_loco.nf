process KINSHIP_LOCO {
    tag "chr${chr}"
    label 'medium'

    publishDir "${params.outdir}/kinship",    mode: 'copy', pattern: "kin_loco_chr*.cXX.txt"
    publishDir "${params.outdir}/chr_bfiles", mode: 'copy', pattern: "chr*.{bed,bim,fam}"

    container params.container

    input:
    tuple val(name), path(bed), path(bim), path(fam), val(chr)
    path pheno

    output:
    tuple val(chr), path("kin_loco_chr${chr}.cXX.txt"),                    emit: kinship
    tuple val(chr), path("chr${chr}.bed"), path("chr${chr}.bim"), path("chr${chr}.fam"), emit: chr_bfiles

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}

    # Chr-excluded bfile → LOCO kinship
    plink --bfile ${name} --not-chr ${chr} --make-bed --out excl_chr${chr} --silent

    gemma -bfile excl_chr${chr} \\
          -p ${pheno} \\
          -n 1 \\
          -gk 1 \\
          -o kin_loco_chr${chr}
    mv output/*.cXX.txt .

    # Chr-only bfile → passed to GWAS step
    plink --bfile ${name} --chr ${chr} --make-bed --out chr${chr} --silent
    """
}
