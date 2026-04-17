# mouse-statgen

A modular Nextflow DSL2 pipeline for GWAS in inbred mouse strains, built around GEMMA, GCTA, and SuSiE.

## Overview

The pipeline takes phenotype data and HMDP genotype data, prepares PLINK files, computes LOCO kinship matrices, runs linear mixed model association tests, fine-maps significant loci with SuSiE, annotates genes, estimates heritability and genetic correlations, and generates publication-quality Manhattan plots.

## Pipeline stages

| # | Stage | Description |
|---|-------|-------------|
| 0 | Prepare inputs | Filter strains/SNPs, QQ-normalize phenotypes, build PLINK and GEMMA files |
| 1 | Kinship | Full-genome and LOCO kinship matrices (GEMMA); LDAK and GCTA GRMs |
| 2 | Heritability | PVE/h² estimation via GEMMA, LDAK, and GCTA |
| 3 | Permutation threshold | Empirical FWER threshold from permuted GWAS |
| 4 | GEMMA LMM | LOCO association testing per trait per chromosome |
| 4b | Combine GWAS | Wide-format effect and p-value matrices across all traits |
| 5 | LD clumping | Per-trait and combined clumping (all traits + per group) |
| 6 | SuSiE fine-mapping | Credible sets and PIPs per locus per trait |
| 7 | Gene annotation | Ensembl biomaRt gene lookup for credible set or lead SNPs |
| 8 | Genetic correlations | Pairwise rg via GCTA bivariate REML |
| 9 | Manhattan plots | Per-trait and per-group plots with credible set labels |
| 10 | Peak clustering | K-means clustering of loci by cross-trait p-value patterns |

## Usage

```bash
NXF_VER=22.04.3 nextflow run main.nf \
    --input      data/anxfb_genetics.csv  \
    --yaml       data/anxfb.yaml          \
    --genotypes  data/muster_hmdp.csv.gz  \
    --name       my_run                   \
    --container  containers/stat-gen.sif  \
    --outdir     output                   \
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
| `--genotypes` | HMDP genotype matrix (gzipped CSV) | — |
| `--name` | Run name prefix for output files | — |
| `--outdir` | Output directory | `results` |
| `--container` | Path to Singularity `.sif` with R/GEMMA/GCTA environment | — |
| `--n_traits` | Number of traits (sets array job size) | `69` |
| `--downsample` | Max individuals per strain (0 = average within strain) | `100000` |
| `--maf` | Min minor allele frequency | `0.1` |
| `--missing` | Max missing genotype rate per SNP | `0.05` |
| `--qqnorm` | QQ-normalise each phenotype before output | `true` |
| `--pval_thresh` | Fallback significance threshold for locus definition | `1e-6` |
| `--perm_trait` | Trait index used for permutation threshold estimation | `1` |
| `--n_perms` | Number of permutations for empirical threshold | `10` |
| `--clump_r2` | LD clumping r² threshold | `0.4` |
| `--clump_kb` | LD clumping window (kb) | `10000` |
| `--susie_L` | Max causal signals per locus in SuSiE | `10` |
| `--locus_window` | Flanking window around lead SNP for fine-mapping (bp) | `1e6` |
| `--annot_input` | SNPs passed to gene annotation: `susie` (credible sets) or `clump` (lead SNPs) | `susie` |
| `--label_gap` | Min distance (kb) between labeled loci on Manhattan plots | `1000` |
| `--k_clusters` | Number of k-means clusters for peak clustering | `7` |
| `--pvalthr_clusters` | −log₁₀(p) threshold for cluster Manhattan plots | `5` |
| `--corr_traits` | Optional file listing trait subset for genetic correlations | — |

## Input files

| File | Description |
|------|-------------|
| Phenotype CSV | Individual-level phenotypes with strain metadata |
| YAML config | Trait definitions (`papername`, `group`), F1 strain mappings, covariate spec, confounding SNPs to exclude |
| Genotype CSV | HMDP genotype matrix (strains × SNPs, gzipped) |
| Correlation traits list | Optional plain-text list of trait names to include in genetic correlation stage |

## Directory structure

```
mouse-statgen/
├── main.nf             # Main workflow
├── nextflow.config     # SLURM + Singularity config
├── submit.sh           # Cluster submission script
├── data/               # Input files (not tracked in git)
├── containers/         # Singularity .sif files (not tracked in git)
├── modules/            # Nextflow process modules
└── src/                # R scripts called by modules
                        # (manhattan.R, manhattan_gg.R, manhattan_local.R,
                        #  postprocess_mv.R, postprocess_pg.R not tracked in git)
```

## References

- Zhou, X. & Stephens, M. Efficient multivariate linear mixed model algorithms for genome-wide association studies. *Nature Methods* **11**, 407–409 (2014). (GEMMA)
- Yang, J., Lee, S.H., Goddard, M.E. & Visscher, P.M. GCTA: a tool for genome-wide complex trait analysis. *The American Journal of Human Genetics* **88**, 76–82 (2011). (GCTA)
- Speed, D. et al. Reevaluation of SNP heritability in complex human traits. *Nature Genetics* **49**, 986–992 (2017). (LDAK)
- Wang, G., Sarkar, A., Carbonetto, P. & Stephens, M. A simple new approach to variable selection in regression, with application to genetic fine mapping. *Journal of the Royal Statistical Society Series B* **82**, 1273–1300 (2020). (SuSiE)
