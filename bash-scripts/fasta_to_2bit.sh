#!/bin/bash

#SBATCH --job-name=fa2bit
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000MB
#SBATCH --ntasks=1

module add SequenceAnalysis/blat/36

# $1 requires the genome reference as FASTA file and $2 .2bit version

faToTwoBit $1 $2
