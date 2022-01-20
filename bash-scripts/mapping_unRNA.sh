#!/bin/bash

#SBATCH --job-name=mapping
#SBATCH --cpus-per-task=4
#SBATCH --mem=8000MB

module add UHTS/Aligner/bowtie/1.2.0;

for x in $(ls -d *tr.fastq); \
do echo ${x}; \
bowtie \
-S \
-t \
-p 4 \
../genome_annotation/tRNA/Rnor_r_sno_sn_t_RNA ${x} \
--un $(basename ${x} .fastq)_no_r_sno_sn_t-RNA.fastq \
2> $(basename ${x} .fastq)_no_r_sno_sn_t-RNA_log.txt > /dev/null; done
