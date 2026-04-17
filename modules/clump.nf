process CLUMP {
    tag "trait${trait_idx}"
    label 'low'

    publishDir "${params.outdir}/clump", mode: 'copy'

    container params.container

    input:
    tuple val(trait_idx), path(assoc)
    tuple val(name), path(bed), path(bim), path(fam)
    path threshold

    output:
    tuple val(trait_idx), path("trait${trait_idx}.clumped"), emit: clumped

    script:
    """
    pthresh=\$(cat ${threshold})

    # Reformat GEMMA output to PLINK clump format (SNP + P columns), using p_lrt
    awk 'NR==1 {for(i=1;i<=NF;i++) if(\$i=="p_lrt") col=i}
         NR>1  {print \$2, \$col}
         NR==1 {print "SNP P"}' ${assoc} > trait${trait_idx}_plink.txt

    plink --bfile ${name} \\
          --clump          trait${trait_idx}_plink.txt \\
          --clump-p1       \$pthresh \\
          --clump-p2       \$pthresh \\
          --clump-r2       ${params.clump_r2} \\
          --clump-kb       ${params.clump_kb} \\
          --out            trait${trait_idx}

    # If no significant SNPs, create empty file so the channel doesn't stall
    [ -f trait${trait_idx}.clumped ] || touch trait${trait_idx}.clumped
    """
}
