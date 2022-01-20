#!/bin/bash

#SBATCH --job-name=sam_summary
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000MB
#SBATCH --nodes=1

module add UHTS/Analysis/samtools/1.10

samtools stats $1
