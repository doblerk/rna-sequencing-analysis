#!/bin/bash

#SBATCH --job-name=fastqc-dump
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module add UHTS/Analysis/sratoolkit/2.10.7;

fastq-dump --gzip SRR*.sra
