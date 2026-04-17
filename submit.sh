#!/bin/bash
#SBATCH -q batch -p compute
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 2:00:00
#SBATCH --mem=4GB
#SBATCH -o slurm.%j.out

# ── Load modules ──────────────────────────────────────────────────────────────
module use --append /projects/omics_share/meta/modules
module load nextflow
module load singularity

cwd=/projects/kumar-lab/USERS/sabnig/

# ── Run pipeline ──────────────────────────────────────────────────────────────
NXF_VER=22.04.3 nextflow run main.nf \
    --input  ${cwd}/pipes/mouse-statgen/data/anxfb_genetics.csv \
    --yaml   ${cwd}/pipes/mouse-statgen/data/anxfb.yaml \
    --name   anxfb_test \
    --outdir output \
    -resume
