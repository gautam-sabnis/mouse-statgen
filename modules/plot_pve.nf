process PLOT_PVE {
    tag "${name}"
    label 'low'

    publishDir "${params.outdir}/plots/heritability", mode: 'copy'

    container params.container

    input:
    path gemma_csv  // heritability_gemma_<name>.csv from COLLECT_HERITABILITY
    path gcta_csv   // heritability_gcta_<name>.csv  from COLLECT_HERITABILITY_GCTA
    path ldak_csv   // heritability_<name>.csv        from HERITABILITY_LDAK
    val  name
    path yaml
    path script     // plot_pve.R
    path aesthetics // aesthetics.R  — staged alongside script so dirname() lookup works
    path setup      // postprocess_setup.R — staged alongside script

    output:
    path "PVE_plot_*.pdf", emit: plots

    script:
    """
    for pair in "gemma ${gemma_csv}" "gcta ${gcta_csv}" "ldak ${ldak_csv}"; do
        method=\$(echo "\$pair" | cut -d' ' -f1)
        csv=\$(   echo "\$pair" | cut -d' ' -f2)
        Rscript ${script} \\
            --input   "\$csv"   \\
            --method  "\$method" \\
            --name    ${name}   \\
            --yaml    ${yaml}   \\
            --outdir  .
    done
    """
}
