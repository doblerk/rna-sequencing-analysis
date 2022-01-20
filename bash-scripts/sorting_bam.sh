#!/bin/bash

#SBATCH --job-name=mapping
#SBATCH --cpus-per-task=4
#SBATCH --mem=8000MB

module add UHTS/Analysis/samtools/1.10;

for x in $(ls -d *.bam); \
do echo ${x}; \
samtools sort \
-@ 4 \
${x} \
-o $(basename ${x} .bam)_sorted.bam; done
