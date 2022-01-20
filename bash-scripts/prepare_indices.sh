#!/bin/bash

#SBATCH --job-name=bowtie
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000MB
#SBATCH --ntasks=1

module add UHTS/Aligner/bowtie/1.2.0;

# $1 requires a FASTA file and $2 requires a indexed file

bowtie-build $1 $2
