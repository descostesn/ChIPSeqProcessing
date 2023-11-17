#!/bin/bash

#SBATCH --nodes 1
#SBATCH --mem=130gb
#SBATCH --job-name testRead
#SBATCH --ntasks=1
#SBATCH --time 04:50:00
#SBATCH --output slurm_%x_%A_%a.out


Rscript bamReadsStats_test.R
