#!/bin/bash

#SBATCH --job-name=trimming
#SBATCH --cpus-per-task=4
#SBATCH --mem=8000MB
#SBATCH --ntasks=1

module add UHTS/Quality_control/cutadapt/2.5;

for x in $(ls -d *_clpd.fastq.gz); do echo ${x}; \
cutadapt \
-j 12 \
-q 25 \
--cut -10 \
--minimum-length 22 \
-o $(basename ${x} .fastq.gz)_tr.fastq.gz \
${x} 1> $(basename ${x} .fastq.gz)_tr_log.txt; done
