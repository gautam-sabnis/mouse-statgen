process HERITABILITY_GCTA {
    tag "trait${trait_idx}"
    label 'low'

    container params.container

    input:
    val  trait_idx
    tuple path(grm_bin), path(grm_N_bin), path(grm_id)
    path pheno_plink
    path covars_plink

    output:
    path "gcta_h2_${trait_idx}.csv", emit: csv

    script:
    """
    gcta --reml \\
         --grm gcta_kinship \\
         --pheno ${pheno_plink} \\
         --mpheno ${trait_idx} \\
         --covar ${covars_plink} \\
         --out gcta_h2 \\
         --thread-num ${task.cpus}

    h2=\$(awk '/^V\\(G\\)\\/Vp/{print \$2}' gcta_h2.hsq)
    se=\$(awk '/^V\\(G\\)\\/Vp/{print \$3}' gcta_h2.hsq)
    printf '%s,%s,%s\\n' "${trait_idx}" "\$h2" "\$se" > gcta_h2_${trait_idx}.csv
    """
}

process COLLECT_HERITABILITY_GCTA {
    tag "${name}"
    label 'low'

    publishDir "${params.outdir}/heritability", mode: 'copy'

    container params.container

    input:
    val  name
    path csvs
    path pheno_order

    output:
    path "heritability_gcta_${name}.csv", emit: heritability

    script:
    """
    mapfile -t trait_names < <(awk 'NR>1 {print \$1}' ${pheno_order})

    printf 'trait_idx,trait,h2,se\\n' > heritability_gcta_${name}.csv

    cat ${csvs} | sort -t',' -k1,1n | while IFS=',' read -r idx h2 se; do
        trait="\${trait_names[\$((idx-1))]}"
        printf '%s,%s,%s,%s\\n' "\$idx" "\$trait" "\$h2" "\$se"
    done >> heritability_gcta_${name}.csv
    """
}
