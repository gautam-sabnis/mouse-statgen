# mouse-statgen

A modular Nextflow DSL2 pipeline for GWAS in inbred mouse strains, built around GEMMA, GCTA, and SuSiE.

## Overview

The pipeline takes phenotype data and HMDP genotype data, prepares PLINK files, runs linear mixed model association tests with LOCO kinship matrices, and fine-maps significant loci.

## Pipeline stages

| # | Module | Status |
|---|--------|--------|
| 0 | Prepare inputs (filter strains/SNPs, build PLINK files) | in progress |
| 1 | GEMMA kinship — LOCO kinship matrices | planned |
| 2 | GEMMA LMM — association testing per trait per chromosome | planned |
| 3 | SuSiE — fine-mapping of significant loci | planned |
| 4 | Manhattan plots — per-trait with credible sets | planned |
| 5 | Group Manhattan — ACAT-combined p-values per phenotype group | planned |
| 6 | Pleiotropy matrix — cross-trait p-values at SuSiE loci | planned |
| 7 | Peak clustering — cluster loci by cross-trait patterns | planned |

## Usage

```bash
NXF_VER=22.04.3 nextflow run main.nf \
    --input      data/anxfb_genetics.csv  \
    --yaml       data/anxfb.yaml          \
    --genotypes  data/muster_hmdp.csv.gz  \
    --name       my_run                   \
    --r_container containers/stat-gen.sif \
    --outdir     output                   \
    -profile slurm,singularity
```

On a SLURM cluster, submit via:

```bash
sbatch submit.sh
```

## Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Phenotype CSV |
| `--yaml` | Run config YAML |
| `--genotypes` | HMDP genotype CSV (gzipped) |
| `--name` | Run name prefix for output files |
| `--outdir` | Output directory (default: `results`) |
| `--r_container` | Path to Singularity `.sif` with R environment |
| `--maf` | Min minor allele frequency (default: 0.1) |
| `--missing` | Max missing genotype rate per SNP (default: 0.05) |
| `--pval_thresh` | Significance threshold for locus definition (default: 1e-6) |

## Directory structure

```
mouse-statgen/
├── main.nf             # Main workflow
├── nextflow.config     # SLURM + Singularity config
├── submit.sh           # Cluster submission script
├── data/               # Input files (not tracked in git)
├── containers/         # Singularity .sif files (not tracked in git)
├── modules/            # Nextflow process modules
└── src/                # R and bash scripts called by modules
```
