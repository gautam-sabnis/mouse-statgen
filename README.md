# mouse-statgen

A modular Nextflow DSL2 pipeline for GWAS in inbred mouse strains, built around GEMMA, GCTA, and SuSiE.

## Overview

The pipeline takes phenotype data and HMDP or MDA genotype data, prepares PLINK files, computes LOCO kinship matrices, runs linear mixed model association tests, fine-maps significant loci (optionally via SuSiE), annotates genes, estimates heritability and genetic correlations, and generates publication-quality Manhattan plots, heritability comparisons, LD decay curves, and INRICH pathway enrichment results.

## Pipeline stages

| # | Stage | Description |
|---|-------|-------------|
| 0 | Prepare inputs | Filter strains/SNPs, QQ-normalize phenotypes, build PLINK and GEMMA files |
| 0b | LD decay plot | Genome-wide LD decay curve from PLINK r² output |
| 1 | Kinship | Full-genome and LOCO kinship matrices (GEMMA); LDAK and GCTA GRMs |
| 2 | Heritability | PVE/h² estimation via GEMMA, LDAK, and GCTA |
| 2c | Heritability plots | Side-by-side PVE comparison plots across GEMMA, GCTA, and LDAK |
| 3 | Permutation threshold | Empirical FWER threshold from permuted GWAS |
| 4 | GEMMA LMM | LOCO association testing per trait per chromosome |
| 4b | Combine GWAS | Wide-format effect and p-value matrices across all traits |
| 5 | LD clumping | Per-trait and combined clumping (all traits + per group) |
| 6 | SuSiE fine-mapping | Credible sets and PIPs per locus per trait (optional; skip with `--run_susie false`) |
| 7 | Gene annotation | Ensembl biomaRt gene lookup for credible set or lead SNPs |
| 8 | Genetic correlations | Pairwise rg via GCTA bivariate REML |
| 9 | Manhattan plots | Per-trait and per-group plots with credible set labels |
| 10 | Peak clustering | K-means clustering of loci by cross-trait p-value patterns |
| 11 | INRICH | Pathway/gene-set enrichment per trait, per group, and per cluster |
| 12 | INRICH plots | Heatmap summaries of INRICH results across gene sets |

## Usage

```bash
NXF_VER=22.04.3 nextflow run main.nf \
    --input      data/phenotypes.csv    \
    --yaml       data/run.yaml          \
    --genotypes  data/muster_hmdp.csv.gz \
    --name       my_run                 \
    --container  containers/stat-gen.sif \
    --outdir     output                 \
    -profile slurm,singularity
```

To use MDA chip genotypes instead of the single HMDP matrix:

```bash
NXF_VER=22.04.3 nextflow run main.nf \
    --input   data/phenotypes.csv  \
    --yaml    data/run.yaml        \
    --mda     true                 \
    --extdata data/extdata         \
    --name    my_run               \
    -profile slurm,singularity
```

On a SLURM cluster, submit via:

```bash
sbatch submit.sh
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Phenotype CSV (individual-level) | — |
| `--yaml` | Run config YAML (trait definitions, groups, covariates) | — |
| `--genotypes` | HMDP genotype matrix (gzipped CSV); ignored when `--mda true` | — |
| `--mda` | Use MDA chip files from `--extdata` instead of a single HMDP matrix | `false` |
| `--extdata` | Directory of per-chromosome MDA genotype files (`snp_*.csv.gz`); required when `--mda true` | — |
| `--name` | Run name prefix for output files | — |
| `--outdir` | Output directory | `results` |
| `--container` | Path to Singularity `.sif` with R/GEMMA/GCTA environment | — |
| `--n_traits` | Number of traits (sets array job size; auto-detected from output if not set) | `69` |
| `--downsample` | Max individuals per strain (0 = average within strain) | `10` |
| `--maf` | Min minor allele frequency | `0.1` |
| `--missing` | Max missing genotype rate per SNP | `0.05` |
| `--qqnorm` | QQ-normalise each phenotype before output | `true` |
| `--pval_thresh` | Fallback significance threshold for locus definition | `1e-6` |
| `--pval_type` | p-value column for association results: `p_lrt`, `p_score`, or `p_wald` | `p_wald` |
| `--perm_trait` | Trait index used for permutation threshold estimation | `1` |
| `--n_perms` | Number of permutations for empirical threshold | `100` |
| `--clump_r2` | LD clumping r² threshold | `0.4` |
| `--clump_kb` | LD clumping window (kb) | `10000` |
| `--run_susie` | Run SuSiE fine-mapping; when `false`, LD clumping is used for all downstream steps | `false` |
| `--susie_L` | Max causal signals per locus in SuSiE | `10` |
| `--locus_window` | Flanking window around lead SNP for fine-mapping (bp) | `1e6` |
| `--loddrop` | LOD-drop threshold for interval expansion around lead SNP (`0` = fixed `locus_window`) | `1.5` |
| `--annot_input` | SNPs passed to gene annotation: `susie` (credible sets) or `clump` (lead SNPs) | `susie` |
| `--label_gap` | Min distance (kb) between labeled loci on Manhattan plots | `1000` |
| `--k_clusters` | Number of k-means clusters for peak clustering | `7` |
| `--pvalthr_clusters` | −log₁₀(p) threshold for cluster Manhattan plots | `5` |
| `--corr_traits` | File listing trait subset for genetic correlations (one name per line) | `data/corr_traits.txt` |
| `--inrich_static` | Directory of pre-generated INRICH gene-set files (from `generate_inrich_static.R`) | `data/inrich_static` |
| `--inrich_i` | Min gene-set size passed to `inrich -i` | `5` |
| `--inrich_j` | Max gene-set size passed to `inrich -j` | `200` |
| `--inrich_minpval` | Adjusted p-value threshold for INRICH heatmap display | `0.05` |

## Input files

| File | Description |
|------|-------------|
| Phenotype CSV | Individual-level phenotypes with strain metadata |
| YAML config | Trait definitions (`papername`, `group`), F1 strain mappings, covariate spec, confounding SNPs to exclude |
| Genotype CSV | HMDP genotype matrix (strains × SNPs, gzipped); or per-chromosome MDA files when `--mda true` |
| Correlation traits list | Plain-text list of trait names (one per line) for the genetic correlation stage |
| INRICH static files | Pre-generated gene-set files in `data/inrich_static/`; generate once with `src/generate_inrich_static.R` |

## INRICH setup

INRICH requires static gene-set and gene-coordinate files that are generated once before running the pipeline:

```bash
Rscript src/generate_inrich_static.R \
    --extdata  data/extdata \
    --outdir   data/inrich_static
```

This requires internet access (biomaRt) and produces files for KEGG pathways, GO terms (BP/CC/MF), mouse phenotype ontology (MP/MGI), brain expression groups, and optionally SFARI gene sets. Commit the outputs to `data/inrich_static/` and point `--inrich_static` at that directory.

## Directory structure

```
mouse-statgen/
├── main.nf             # Main workflow
├── nextflow.config     # SLURM + Singularity config
├── submit.sh           # Cluster submission script
├── data/               # Input files (not tracked in git)
│   └── inrich_static/  # Pre-generated INRICH gene-set files
├── containers/         # Singularity .sif files (not tracked in git)
├── modules/            # Nextflow process modules
│   ├── prepare_input.nf
│   ├── kinship_{full,loco,ldak,gcta}.nf
│   ├── heritability{,_ldak,_gcta}.nf
│   ├── gemma_lmm.nf
│   ├── permutation.nf
│   ├── clump.nf / clump_combined.nf
│   ├── susie_finemap.nf
│   ├── get_genes.nf
│   ├── combine_gwas.nf
│   ├── genetic_corr.nf
│   ├── plot_manhattan.nf
│   ├── cluster_peaks.nf
│   ├── plot_pve.nf
│   ├── plot_ld_decay.nf
│   ├── inrich.nf
│   └── plot_inrich.nf
└── src/                # R scripts called by modules
    ├── prepare_and_make_plink.R
    ├── susie_finemap.R
    ├── plot_manhattan.R
    ├── aesthetics.R
    ├── postprocess_setup.R
    ├── cluster_peaks.R
    ├── plot_pve.R
    ├── plot_ld_decay.R
    ├── write_inrich_intervals.R
    ├── write_inrich_intervals_combined.R
    ├── generate_inrich_static.R  # run once, offline
    └── plot_inrich.R
```

## References

- Zhou, X. & Stephens, M. Efficient multivariate linear mixed model algorithms for genome-wide association studies. *Nature Methods* **11**, 407–409 (2014). (GEMMA)
- Yang, J., Lee, S.H., Goddard, M.E. & Visscher, P.M. GCTA: a tool for genome-wide complex trait analysis. *The American Journal of Human Genetics* **88**, 76–82 (2011). (GCTA)
- Speed, D. et al. Reevaluation of SNP heritability in complex human traits. *Nature Genetics* **49**, 986–992 (2017). (LDAK)
- Wang, G., Sarkar, A., Carbonetto, P. & Stephens, M. A simple new approach to variable selection in regression, with application to genetic fine mapping. *Journal of the Royal Statistical Society Series B* **82**, 1273–1300 (2020). (SuSiE)
- Lee, I. et al. A single nucleotide polymorphism tagging set for genetic association studies in the mouse. *Genomics* **89**, 104–107 (2007). (INRICH)
