process COMBINE_GWAS {
    tag "${name}"
    label 'medium'

    publishDir "${params.outdir}/gwas", mode: 'copy'

    container params.container

    input:
    val  name
    path assoc_files, stageAs: "assoc/*"
    path pheno_order

    output:
    path "gwas_effects_${name}.csv", emit: effects
    path "gwas_pvals_${name}.csv",   emit: pvals

    script:
    """
    # ── Trait names (1-based, skip header) ───────────────────────────────────
    mapfile -t trait_names < <(awk 'NR>1 {print \$1}' ${pheno_order})

    # ── Sort assoc files by trait index ──────────────────────────────────────
    mapfile -t files < <(ls assoc/trait*_loco.assoc.txt | sort -V)

    # ── Find column indices from header ───────────────────────────────────────
    header=\$(head -1 "\${files[0]}")
    chr_col=\$(echo "\$header"  | tr '\\t' '\\n' | grep -n '^chr\$'     | cut -d: -f1)
    rs_col=\$(echo "\$header"   | tr '\\t' '\\n' | grep -n '^rs\$'      | cut -d: -f1)
    ps_col=\$(echo "\$header"   | tr '\\t' '\\n' | grep -n '^ps\$'      | cut -d: -f1)
    a1_col=\$(echo "\$header"   | tr '\\t' '\\n' | grep -n '^allele1\$' | cut -d: -f1)
    a0_col=\$(echo "\$header"   | tr '\\t' '\\n' | grep -n '^allele0\$' | cut -d: -f1)
    beta_col=\$(echo "\$header" | tr '\\t' '\\n' | grep -n '^beta\$'    | cut -d: -f1)
    se_col=\$(echo "\$header"   | tr '\\t' '\\n' | grep -n '^se\$'      | cut -d: -f1)
    plrt_col=\$(echo "\$header" | tr '\\t' '\\n' | grep -n '^p_lrt\$'   | cut -d: -f1)
    psco_col=\$(echo "\$header" | tr '\\t' '\\n' | grep -n '^p_score\$' | cut -d: -f1)

    # ── Extract SNP key columns (space-separated); recode chr 23 → X ────────────
    awk -v c="\$chr_col" -v r="\$rs_col" -v p="\$ps_col" \\
        -v a1="\$a1_col" -v a0="\$a0_col" \\
        'NR>1 {OFS=" "; chr=(\$c=="23"?"X":\$c); print chr,\$r,\$p,\$a1,\$a0}' "\${files[0]}" > keys.tmp

    # ── Loop over trait files ─────────────────────────────────────────────────
    eff_files=(keys.tmp)
    pval_files=(keys.tmp)
    eff_hdr="chr rs ps allele1 allele0"
    pval_hdr="chr rs ps allele1 allele0"

    for i in "\${!files[@]}"; do
        f="\${files[\$i]}"
        tidx=\$(basename "\$f" | sed 's/trait\\([0-9]*\\)_loco\\.assoc\\.txt/\\1/')
        tname="\${trait_names[\$((tidx-1))]}"

        awk -v b="\$beta_col" -v s="\$se_col" \\
            'NR>1 {OFS=" "; print \$b,\$s}' "\$f" > eff_\${i}.tmp
        awk -v pl="\$plrt_col" -v ps="\$psco_col" \\
            'NR>1 {OFS=" "; print \$pl,\$ps}' "\$f" > pval_\${i}.tmp

        eff_files+=(eff_\${i}.tmp)
        pval_files+=(pval_\${i}.tmp)
        eff_hdr="\${eff_hdr} beta_\${tname} se_\${tname}"
        pval_hdr="\${pval_hdr} p_lrt_\${tname} p_score_\${tname}"
    done

    # ── Write header + pasted data ────────────────────────────────────────────
    echo "\$eff_hdr"               >  gwas_effects_${name}.csv
    paste -d' ' "\${eff_files[@]}" >> gwas_effects_${name}.csv

    echo "\$pval_hdr"               >  gwas_pvals_${name}.csv
    paste -d' ' "\${pval_files[@]}" >> gwas_pvals_${name}.csv
    """
}
