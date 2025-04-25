#!/bin/bash
#SBATCH --job-name=render_rmd
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --partition=caslake
#SBATCH --account=pi-xinhe

module load pandoc

Rscript render_rmd.R
