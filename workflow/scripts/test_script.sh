#!/bin/bash

#SBATCH --nodes 1
#SBATCH --mem=80gb
#SBATCH --job-name testRead
#SBATCH --ntasks=1
#SBATCH --time 04:50:00
#SBATCH --output slurm_%x_%A_%a.out


conda activate R-4.3.2

Rscript bamReadsStats_test.R

conda deactivate
