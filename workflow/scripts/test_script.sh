#!/bin/bash

#SBATCH --nodes 1
#SBATCH --mem=100gb
#SBATCH --job-name testRead
#SBATCH --ntasks=1
#SBATCH --time 40:50:00
#SBATCH --cpus-per-task=8
#SBATCH --output slurm_%x_%A_%a.out


Rscript bamReadsStats_test_reducebyyield_severalfiles.R
