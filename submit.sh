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
YAML=${cwd}pipes/mouse-statgen/data/grooming_nowild.yaml
N_TRAITS=$(awk '/^phenotypes:/{flag=1; next} flag && /^  - /{count++} flag && /^[^ ]/{flag=0} END{print count}' "${YAML}")


# ── Run pipeline ──────────────────────────────────────────────────────────────
NXF_VER=22.04.3 nextflow run main.nf \
    --input  ${cwd}pipes/mouse-statgen/data/grooming_paper_strain_survey_2019_11_21.csv \
    --yaml   ${YAML} \
    --n_traits ${N_TRAITS} \
    --name   grooming_paper \
    --outdir output_grooming \
    -resume
