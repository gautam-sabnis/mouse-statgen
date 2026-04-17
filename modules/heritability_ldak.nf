process HERITABILITY_LDAK {
    tag "${name}"
    label 'medium'

    publishDir "${params.outdir}/heritability", mode: 'copy'

    container params.container

    input:
    val  name
    path pheno_plink
    path covars_plink
    tuple path(grm_bin), path(grm_id), path(grm_details)
    path pheno_order

    output:
    path "heritability_${name}.csv", emit: heritability

    script:
    """
    ldak --reml ldak_h2 \\
         --pheno ${pheno_plink} \\
         --grm ldak_kinship \\
         --mpheno ALL \\
         --dentist YES \\
         --covar ${covars_plink} \\
         --kinship-details NO

    # Extract trait names from phenotypes_order.txt (single-column, skip header)
    mapfile -t trait_names < <(awk 'NR>1 {print \$1}' ${pheno_order})

    printf 'trait_idx,trait,h2,se\\n' > heritability_${name}.csv

    n=\$(ls ldak_h2.*.reml | wc -l)
    for i in \$(seq 1 \$n); do
        h2=\$(awk '/^Her_All/{print \$2}' ldak_h2.\$i.reml)
        se=\$(awk '/^Her_All/{print \$3}' ldak_h2.\$i.reml)
        trait="\${trait_names[\$((i-1))]}"
        printf '%s,%s,%s,%s\\n' "\$i" "\$trait" "\$h2" "\$se"
    done >> heritability_${name}.csv
    """
}
