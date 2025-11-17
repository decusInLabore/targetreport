#!/bin/bash
#SBATCH --job-name=MAR
#SBATCH --output=quarto_render.out
#SBATCH --error=quarto_render.err
#SBATCH --time=06:00:00
#SBATCH --mem=200G

module load Singularity/3.6.4

singularity shell --cleanenv --bind /nemo:/nemo,/camp:/camp /nemo/stp/babs/working/boeings/projects/boeings/stefan.boeing/619_demo_deepmap/envs/r450.python310.ubuntu.22.04_latest.sif bash -c "
    source /nemo/stp/babs/working/boeings/projects/boeings/stefan.boeing/619_demo_deepmap/envs/demo_venv_310/bin/activate
    Rscript quorto.R
"


