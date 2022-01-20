#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000MB
#SBATCH --ntasks=1

module add UHTS/Quality_control/fastqc/0.11.9;

fastqc -o ./fastqc_clpd_tr  --extract *clpd_tr.fastq.gz --threads 1
