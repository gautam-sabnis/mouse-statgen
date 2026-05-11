// ─────────────────────────────────────────────────────────────────────────────
// INRICH enrichment analysis modules
//
// Four processes:
//   WRITE_INRICH_SNP_MAP          — BIM → SNPs_map_for_INRICH.txt (once)
//   WRITE_INRICH_INTERVALS        — per-trait interval files (parallelised)
//   WRITE_INRICH_INTERVALS_COMBINED — per-group + per-cluster (once)
//   RUN_INRICH                    — runs inrich for one label × all gene-sets
// ─────────────────────────────────────────────────────────────────────────────

// ── SNP map ──────────────────────────────────────────────────────────────────
process WRITE_INRICH_SNP_MAP {
    tag "snp_map"
    label 'low'

    publishDir "${params.outdir}/inrich/inputs", mode: 'copy'

    container params.container

    input:
    tuple val(name), path(bed), path(bim), path(fam)

    output:
    path "SNPs_map_for_INRICH.txt", emit: snp_map

    script:
    """
    awk '{print \$1, \$4}' ${bim} > SNPs_map_for_INRICH.txt
    """
}

// ── Per-trait intervals ───────────────────────────────────────────────────────
process WRITE_INRICH_INTERVALS {
    tag "trait${trait_idx}"
    label 'low'

    container params.container

    input:
    tuple val(trait_idx), path(locus_files)
    path  pheno_order
    val   mode   // "susie" or "clump"
    tuple val(plink_name), path(bed), path(bim), path(fam)

    output:
    path "intervals_*_for_INRICH.txt", emit: intervals, optional: true

    script:
    """
    Rscript ${projectDir}/src/write_inrich_intervals.R \\
        --trait_idx   ${trait_idx}                     \\
        --pheno_order ${pheno_order}                   \\
        --mode        ${mode}                          \\
        --loddrop     ${params.loddrop}                \\
        --maxdist     ${params.clump_kb * 1000}        \\
        --pval_type   ${params.pval_type}              \\
        --plink_prefix ${plink_name}                   \\
        --clump_kb    ${params.clump_kb}               \\
        --clump_r2    ${params.clump_r2}
    """
}

// ── Per-group + per-cluster intervals ─────────────────────────────────────────
process WRITE_INRICH_INTERVALS_COMBINED {
    tag "combined"
    label 'low'

    container params.container

    input:
    path  locus_files     // all trait*_locus_*.txt (SuSiE) or trait*.clumped (clump mode)
    path  pheno_order
    path  trait_groups    // trait_groups_<name>.tsv from PREPARE_INPUT
    path  cluster_files   // *.csv from CLUSTER_PEAKS (includes cluster_assignments_*.csv)
    val   mode            // "susie" or "clump"
    path  assoc_files     // trait*_loco.assoc.txt from COLLECT_GEMMA_LMM (all traits, collected)
    tuple val(plink_name), path(bed), path(bim), path(fam)

    output:
    path "intervals_*_for_INRICH.txt", emit: intervals, optional: true

    script:
    """
    Rscript ${projectDir}/src/write_inrich_intervals_combined.R \\
        --pheno_order  ${pheno_order}              \\
        --trait_groups ${trait_groups}             \\
        --name         ${params.name}              \\
        --mode         ${mode}                     \\
        --loddrop      ${params.loddrop}           \\
        --maxdist      ${params.clump_kb * 1000}   \\
        --pval_type    ${params.pval_type}         \\
        --plink_prefix ${plink_name}               \\
        --clump_kb     ${params.clump_kb}          \\
        --clump_r2     ${params.clump_r2}
    """
}

// ── Run INRICH ────────────────────────────────────────────────────────────────
process RUN_INRICH {
    tag "${label}"
    label 'high'

    publishDir "${params.outdir}/inrich/results", mode: 'copy'

    container params.container

    input:
    tuple val(label), path(interval_file)
    path  snp_map
    path  static_files   // genes_coordinates + all groups_* files from inrich_static/

    output:
    path "*.out.inrich", emit: results, optional: true

    script:
    """
    # Skip if no intervals were written for this label
    if [ ! -s "${interval_file}" ]; then
        echo "Empty interval file for '${label}' — skipping INRICH."
        exit 0
    fi

    for geneset_file in groups_*_for_INRICH.txt; do
        [ -f "\$geneset_file" ] || continue
        geneset=\$(echo "\$geneset_file" | sed 's/^groups_//; s/_for_INRICH\\.txt\$//')
        out_stem="intervals_${label}_groups_\${geneset}"

        inrich -c \\
            -a "${interval_file}"               \\
            -m "${snp_map}"                     \\
            -g genes_coordinates_for_INRICH.txt \\
            -t "\$geneset_file"                 \\
            -o "\$out_stem"                     \\
            -i ${params.inrich_i}               \\
            -j ${params.inrich_j} || true
    done
    """
}
