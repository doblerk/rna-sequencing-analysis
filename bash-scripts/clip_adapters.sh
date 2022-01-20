#!/bin/bash

#SBATCH --job-name=clipping
#SBATCH --cpus-per-task=4
#SBATCH --mem=8000MB
#SBATCH --ntasks=1

module add UHTS/Quality_control/cutadapt/2.5;

for x in $(ls -d *fastq.gz); do echo ${x}; \
cutadapt \
-j 12 \
-a AGATCGGAAGAGCACACGTCTGAA \
-q 25 \
--cut 2 \
--minimum-length 22 \
--discard-untrimmed \
--overlap 3 \
-e 0.2 \
-o $(basename ${x} .fastq.gz)_clpd.fastq.gz \
${x} 1> $(basename ${x} .fastq.gz)_clpd_log.txt; done
