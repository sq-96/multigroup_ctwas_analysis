#!/bin/bash
#SBATCH --job-name=render_rmd
#SBATCH --mem=30G
#SBATCH --cpus-per-task=10
#SBATCH --partition=caslake
#SBATCH --account=pi-guiming

module load pandoc

Rscript render_rmd.R
