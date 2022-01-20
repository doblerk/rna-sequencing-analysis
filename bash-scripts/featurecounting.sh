#!/bin/bash

#SBATCH --job-name=featurecounts
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G

module add UHTS/Analysis/subread/2.0.1

# If you want to counts the reads on CDS:
# -t CDS -g gene_id
featureCounts -T 4 -t exon -g gene_biotype -a ../genome_annotation/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.104.gtf -o biotype_counts_rawfile.txt *_sorted.bam 
