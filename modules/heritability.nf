process HERITABILITY {
    tag "trait${trait_idx}"
    label 'low'

    container params.container

    input:
    tuple val(name), path(bed), path(bim), path(fam), val(trait_idx)
    path pheno
    path covars
    path kinship

    output:
    path "heritability_${trait_idx}.tsv", emit: tsv

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}

    awk -F'\\t' -v col=${trait_idx} '{print (\$col=="" ? "NA" : \$col)}' ${pheno} > single_pheno.txt

    gemma -bfile ${name} \\
          -p single_pheno.txt \\
          -n 1 \\
          -k ${kinship} \\
          -c ${covars} \\
          -lmm \\
          -o run

    log=output/run.log.txt
    pve=\$(  grep "^## pve estimate"  \$log | awk '{print \$NF}')
    pvese=\$(grep "^## se(pve)"       \$log | awk '{print \$NF}')
    vg=\$(   grep "^## vg estimate"   \$log | awk '{print \$NF}')
    ve=\$(   grep "^## ve estimate"   \$log | awk '{print \$NF}')

    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "${trait_idx}" "\$pve" "\$pvese" "\$vg" "\$ve" \\
        > heritability_${trait_idx}.tsv
    """
}

process COLLECT_HERITABILITY {
    tag "${name}"
    label 'low'

    publishDir "${params.outdir}/heritability", mode: 'copy'

    container params.container

    input:
    val  name
    path tsvs
    path pheno_order

    output:
    path "heritability_gemma_${name}.csv", emit: heritability

    script:
    """
    mapfile -t trait_names < <(awk 'NR>1 {print \$1}' ${pheno_order})

    printf 'trait_idx,trait,PVE,PVESE,Vg,Ve\\n' > heritability_gemma_${name}.csv

    cat ${tsvs} | sort -k1,1n | while IFS=\$'\\t' read -r idx pve pvese vg ve; do
        trait="\${trait_names[\$((idx-1))]}"
        printf '%s,%s,%s,%s,%s,%s\\n' "\$idx" "\$trait" "\$pve" "\$pvese" "\$vg" "\$ve"
    done >> heritability_gemma_${name}.csv
    """
}
