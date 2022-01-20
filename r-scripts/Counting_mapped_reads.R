#---------------------------------------
# Counting mapped reads
#
#
# This script comprises hands-on code
#---------------------------------------

library(Rsubread)

genome_bam <- c("Biever_Neuropil_Poly_1_clpd_tr_no_r_sno_sn_t-RNA_Rnor_6_sorted.bam",
                "Biever_Neuropil_Poly_2_clpd_tr_no_r_sno_sn_t-RNA_Rnor_6_sorted.bam",
                "Biever_Neuropil_Poly_3_clpd_tr_no_r_sno_sn_t-RNA_Rnor_6_sorted.bam",
                "Biever_Somata_Poly_1_clpd_tr_no_r_sno_sn_t-RNA_Rnor_6_sorted.bam",
                "Biever_Somata_Poly_2_clpd_tr_no_r_sno_sn_t-RNA_Rnor_6_sorted.bam",
                "Biever_Somata_Poly_3_clpd_tr_no_r_sno_sn_t-RNA_Rnor_6_sorted.bam")

genome_annot <- "Rattus_norvegicus.Rnor_6.0.104.gtf"

# Retrieving CDS
fc.matrix <- featureCounts(files = genome_bam,
                           annot.inbuilt = "mm10",
                           annot.ext = genome_annot,
                           isGTFAnnotationFile = TRUE,
                           GTF.featureType = "CDS",
                           GTF.attrType = "gene_id",
                           )

# Retrieving biotypes
fc.matrix <- featureCounts(files = genome_bam,
                           annot.inbuilt = "mm10",
                           annot.ext = genome_annot,
                           isGTFAnnotationFile = TRUE,
                           GTF.featureType = "exon",
                           GTF.attrType = "gene_biotype",
                           )

# Saving the reads count to a .txt file
write.table(fc.matrix$counts, "name of the file", sep = " ", row.names = TRUE, col.names = TRUE)