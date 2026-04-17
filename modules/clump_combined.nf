process CLUMP_COMBINED {
    tag "${name}"
    label 'medium'

    publishDir "${params.outdir}/clump_combined", mode: 'copy'

    container params.container

    input:
    path  pvals          // gwas_pvals_<name>.csv  (space-separated) from COMBINE_GWAS
    path  trait_groups   // trait_groups_<name>.tsv (trait_name\tgroup) from PREPARE_INPUT
    tuple val(name), path(bed), path(bim), path(fam)
    path  threshold      // single-number file from CALC_PERM_THRESHOLD

    output:
    path "*.clumped", emit: clumped   // all.clumped + one per group

    script:
    """
    pthresh=\$(cat ${threshold})

    # ── All-traits combined: row-wise min over every p_lrt_* column ──────────
    awk 'NR==1 {
             for (i=1; i<=NF; i++) {
                 if (\$i == "rs")        rs_col = i
                 if (\$i ~ /^p_lrt_/)   pcols[i] = 1
             }
             print "SNP P"
         }
         NR>1 {
             min_p = 1
             for (i in pcols)
                 if (\$i != "NA" && \$i+0 < min_p) min_p = \$i+0
             if (min_p < 1) print \$rs_col, min_p
         }' ${pvals} > all.plink.txt

    plink --bfile ${name} \\
          --clump       all.plink.txt \\
          --clump-p1    \$pthresh \\
          --clump-p2    \$pthresh \\
          --clump-r2    ${params.clump_r2} \\
          --clump-kb    ${params.clump_kb} \\
          --out         all
    [ -f all.clumped ] || touch all.clumped

    # ── Per-group combined: row-wise min over that group's p_lrt_* columns ───
    # Skip TSV header; each remaining line is "<trait_name>\t<group>"
    tail -n +2 ${trait_groups} | cut -f2 | sort -u | while IFS= read -r grp; do

        # Filesystem-safe stem: replace non-alphanumeric runs with '_', trim trailing '_'
        safe=\$(printf '%s' "\$grp" | tr -cs 'A-Za-z0-9' '_' | sed 's/_\$//')

        # Build PLINK input: SNP + min p_lrt across traits belonging to this group
        awk -v grp="\$grp" -v tsvfile="${trait_groups}" '
            BEGIN {
                while ((getline line < tsvfile) > 0) {
                    n = split(line, a, "\t")
                    if (n == 2 && a[2] == grp) traits[a[1]] = 1
                }
            }
            NR==1 {
                for (i=1; i<=NF; i++) {
                    if (\$i == "rs") rs_col = i
                    col = \$i
                    sub(/^p_lrt_/, "", col)
                    if (col in traits) pcols[i] = 1
                }
                print "SNP P"
            }
            NR>1 {
                min_p = 1
                for (i in pcols)
                    if (\$i != "NA" && \$i+0 < min_p) min_p = \$i+0
                if (min_p < 1) print \$rs_col, min_p
            }
        ' ${pvals} > "\${safe}.plink.txt"

        # Skip groups that produced no SNPs (e.g. all p-values were NA)
        [ -s "\${safe}.plink.txt" ] || { touch "\${safe}.clumped"; continue; }

        plink --bfile ${name} \\
              --clump       "\${safe}.plink.txt" \\
              --clump-p1    \$pthresh \\
              --clump-p2    \$pthresh \\
              --clump-r2    ${params.clump_r2} \\
              --clump-kb    ${params.clump_kb} \\
              --out         "\$safe"
        [ -f "\${safe}.clumped" ] || touch "\${safe}.clumped"

    done
    """
}
