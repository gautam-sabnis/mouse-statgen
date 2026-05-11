process PLOT_LD_DECAY {
    tag "${name}"
    label 'medium'

    publishDir "${params.outdir}/plots", mode: 'copy'

    container params.container

    input:
    tuple val(name), path(bed), path(bim), path(fam)
    path  script     // plot_ld_decay.R
    path  aesthetics // aesthetics.R  — staged alongside script so dirname() lookup works
    path  setup      // postprocess_setup.R — staged alongside script

    output:
    path "LD_decay_${name}.pdf", emit: plot

    script:
    """
    plink --bfile        ${name}  \\
          --chr           1-19    \\
          --r2                    \\
          --ld-window-kb  2500    \\
          --ld-window     999999  \\
          --ld-window-r2  0       \\
          --out           ld_all

    Rscript ${script} \\
        --ld     ld_all.ld \\
        --name   ${name}   \\
        --outdir .
    """
}
