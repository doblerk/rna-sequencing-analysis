#---------------------------------------
# Analysis of biotypes' origin
#
#
# This script comprises hands-on code
#---------------------------------------

library(ggplot2)
library(reshape2)

biotypes <- read.table("biotype_counts_processed.txt", header = T)
colnames(biotypes) <- c("Geneid",paste0(rep(c("Neuropil_Poly_", "Somata_Poly_"), each = 3),1:3))
rownames(biotypes) <- biotypes$Geneid
head(biotypes)

keep <- rowSums(biotypes[,2:ncol(biotypes)]) > 0
biotypes <- biotypes[keep,]

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

biotypes_normalized <- as.data.frame(lapply(biotypes[2:ncol(biotypes)], min_max_norm))
biotypes_normalized$Geneid <- biotypes$Geneid

biotypes_reshaped <- melt(biotypes_normalized, id.vars = c("Geneid"))

pdf("RPF_biotypes.pdf")
ggplot(data = biotypes_reshaped, aes(fill = Geneid, x = variable, y = value)) +
  geom_bar(position="fill", stat="identity") +
  ylab("Proportion of mapped reads") +
  xlab("") +
  scale_fill_discrete(name = "Biotypes") +
  theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.title=element_text(size=12),
        legend.position = "right")
dev.off()