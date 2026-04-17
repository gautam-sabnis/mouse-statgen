process PERMUTATION {
    tag "perm${perm_idx}"
    label 'medium'

    container params.container

    input:
    tuple val(name), path(bed), path(bim), path(fam), val(perm_idx)
    path pheno
    path covars
    path kinship
    val  trait_idx

    output:
    path "minp_${perm_idx}.txt", emit: minp

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}

    # Extract and shuffle the target trait column
    awk -F'\\t' -v col=${trait_idx} 'BEGIN{srand(${perm_idx})} {print (\$col=="" ? "NA" : \$col)}' ${pheno} \\
        | awk 'BEGIN{srand(${perm_idx})} {lines[NR]=\$0} END{
            n=NR; for(i=n;i>1;i--){j=int(rand()*i)+1; tmp=lines[i]; lines[i]=lines[j]; lines[j]=tmp}
            for(i=1;i<=n;i++) print lines[i]
        }' > shuffled_pheno.txt

    gemma -bfile ${name} \\
          -p shuffled_pheno.txt \\
          -n 1 \\
          -k ${kinship} \\
          -c ${covars} \\
          -lmm 4 \\
          -lmin 0.01 -lmax 100 \\
          -o perm${perm_idx}

    # Extract genome-wide minimum p_lrt from this permutation
    awk 'NR==1 {for(i=1;i<=NF;i++) if(\$i=="p_lrt") col=i}
         NR>1  {if(\$col < min || NR==2) min=\$col}
         END   {print min}' output/perm${perm_idx}.assoc.txt \\
        > minp_${perm_idx}.txt
    """
}

process CALC_PERM_THRESHOLD {
    tag "threshold"
    label 'low'

    publishDir "${params.outdir}/gwas", mode: 'copy'

    container params.container

    input:
    path minp_files

    output:
    path "perm_threshold.txt", emit: threshold

    script:
    """
    # Collect all minP values, sort, take 5th percentile (FWER = 0.05)
    cat ${minp_files} | sort -g > all_minp.txt
    n=\$(wc -l < all_minp.txt)
    idx=\$(awk -v n=\$n 'BEGIN{printf "%d", int(n * 0.05) + 1}')
    thresh=\$(sed -n "\${idx}p" all_minp.txt)
    echo "\$thresh" > perm_threshold.txt
    echo "Permutation threshold (FWER=0.05, n=${params.n_perms} perms): \$thresh"
    """
}
