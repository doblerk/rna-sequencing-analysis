#!/bin/bash

#SBATCH --job-name=loading_samples
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G

module add UHTS/Analysis/sratoolkit/2.10.7;

# samples.txt is a txt file containing the BioProject accession number (SRA...)

prefetch --option-file samples.txt
