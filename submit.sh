#!/bin/bash
#SBATCH -q batch -p compute
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 2:00:00
#SBATCH --mem=8GB
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

# ── Load modules ──────────────────────────────────────────────────────────────
module use --append /projects/omics_share/meta/modules
module load nextflow
module load singularity

# ── Run pipeline ──────────────────────────────────────────────────────────────
NXF_VER=22.04.3 nextflow run main.nf \
    --input      data/anxfb_genetics.csv       \
    --yaml       data/anxfb.yaml              \
    --genotypes  data/muster_hmdp.csv.gz  \
    --name       anxfb_test                   \
    --r_container containers/stat-gen.sif      \
    --outdir     output                       \
    -profile slurm,singularity                \
    -resume
