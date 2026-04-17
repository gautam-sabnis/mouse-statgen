process GEMMA_LMM {
    tag "trait${trait_idx}_chr${chr}"
    label 'medium'

    container params.container

    input:
    tuple val(trait_idx), val(chr),
          path(bed), path(bim), path(fam),
          path(kinship)
    path pheno
    path covars

    output:
    tuple val(trait_idx), path("trait${trait_idx}_chr${chr}.assoc.txt"), emit: assoc

    script:
    def bfile = bed.baseName
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}

    awk -F'\\t' -v col=${trait_idx} '{print (\$col=="" ? "NA" : \$col)}' ${pheno} > single_pheno.txt

    gemma -bfile ${bfile} \\
          -k ${kinship} \\
          -lmm 4 \\
          -lmin 0.01 -lmax 100 \\
          -p single_pheno.txt \\
          -n 1 \\
          -c ${covars} \\
          -o trait${trait_idx}_chr${chr}

    mv output/trait${trait_idx}_chr${chr}.assoc.txt .
    """
}

process COLLECT_GEMMA_LMM {
    tag "trait${trait_idx}"
    label 'low'

    publishDir "${params.outdir}/gwas", mode: 'copy'

    container params.container

    input:
    tuple val(trait_idx), path(assoc_files)

    output:
    tuple val(trait_idx), path("trait${trait_idx}_loco.assoc.txt"), emit: assoc

    script:
    """
    # Sort files by chromosome number, write header once then all data rows
    files=\$(ls trait${trait_idx}_chr*.assoc.txt | sort -V)
    head -1 \$(echo \$files | awk '{print \$1}') > trait${trait_idx}_loco.assoc.txt
    for f in \$files; do
        tail -n +2 \$f >> trait${trait_idx}_loco.assoc.txt
    done
    """
}
