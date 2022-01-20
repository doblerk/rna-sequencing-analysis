#---------------------------------------
# Differential expression analysis
#
#
# This script comprises hands-on code
#---------------------------------------

library(DESeq2)
library(gplots)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(Glimma)
library(reshape2)
library(pheatmap)
library(genefilter)
library(gridExtra)
library(org.Rn.eg.db)
library(topGO)


#---------------------------------------
# Constructing a DESeqDataSeq object
#---------------------------------------

# Loading the data
cts <- read.table("CDS_counts.txt", header = T)

# Renaming the data
sample.names <- c(rep(paste0("Neuropil_Poly_",1:3)),
                  rep(paste0("Somata_Poly_",1:3)))

colnames(cts) <- sample.names

# Creating the information about the samples
coldata <- data.frame(c(rep("Neuropil_Poly",3),rep("Somata_Poly",3)),
                      c(rep(c("A","B","C"), 2)))

colnames(coldata) <- c("sample", "comparison")
rownames(coldata) <- c(rep(paste0("Neuropil_Poly_",1:3)),rep(paste0("Somata_Poly_",1:3)))

# Constructing a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ sample)


#---------------------------------------
# Pre-filtering and processing
#---------------------------------------

# Filtering out genes with low read counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Factor levels
dds$condition <- factor(c(rep("Neuropil_Poly",3),rep("Somata_Poly",3)),
                        levels = c("Neuropil_Poly", "Somata_Poly"))


#---------------------------------------
# Differential expression analysis
#---------------------------------------

# Performing differential expression analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("sample", "Neuropil_Poly", "Somata_Poly"), alpha = 0.05)


#---------------------------------------
# Exploratory data analysis
#---------------------------------------

rld <- rlog(dds, blind = F)

# Sample distance matrix
sample.dists <- dist(t(assay(rld)), method = "euclidean")
sample.distsMatrix <- as.matrix(sample.dists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
par(mar=c(0,0,0,5))
pdf("sample_dist_matrix.pdf")
pheatmap(sample.distsMatrix,
         clustering_distance_rows=sample.dists,
         clustering_distance_cols=sample.dists,
         col=colors,
         fontsize_row = 10,
         fontsize_col = 10)
dev.off()

# Visualisation of expression values
rld_reshaped <- melt(assay(rld))
ggplot(data = rld_reshaped, aes(x = Var2, y = value)) +
  geom_violin(aes(fill = Var2) ,trim = T) +
  ggtitle("Violin plots of the expression values") +
  theme(axis.text.x = element_text(angle = 90, size = 12),
        legend.position = 'none', plot.title = element_text(size =12))

# PCA
pca <- prcomp(t(assay(rld)), scale = F)
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
pca$x[,c(1,2)]
pca_df <- as.data.frame(pca$x[,c(1,2)])
pca_df$group <- c(rep("Neuropil_Poly",3),rep("Somata_Poly",3))
par(mar=c(0,0,0,5))
pdf("pca.pdf")
ggplot(data = pca_df, aes(x=PC1, y=PC2)) + 
  geom_point(aes(x=PC1,y=PC2, colour=factor(group)), size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       color = "Groups") +
  scale_fill_discrete(name = "Groups") +
  scale_color_manual(values = c("darkorchid", "royalblue")) +
  theme(axis.text.x = element_text(hjust = 1),
        axis.title=element_text(size=12))
dev.off()

# Expression pheatmap
heat_colors <- brewer.pal(6, "YlOrRd")
annot <- as.data.frame(colData(dds)[, c("condition", "comparison")])
select <- order(rowMeans(counts(dds,normalized = TRUE)))
pdf("pheatmap.pdf")
pheatmap(assay(rld)[select,], cluster_rows = F, cluster_cols = T, show_rownames = F,
         fontsize_col = 12, annotation = annot, fontsize_row = 12)
dev.off()


#---------------------------------------
# Diagnostic plots
#---------------------------------------

# MA plot
plotMA(res, ylim = c(-3, 3))

# Histogram of p-values
hist(res$pvalue[res$baseMean > 1], breaks = 0:30/30,
     col = "grey50", border = "white", xlab = "",
     main = "Histogram of p values for genes with mean normalized count > 1")

# Dispersion plot
plotDispEsts(dds, ylim = c(1e-6, 1e1))

# Gene clustering
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 30)
mat  <- assay(rld)[topVarGenes,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(dds)[, c("condition", "comparison")])
pheatmap(mat, annotation_col = anno)


#---------------------------------------
# Gene annotation
#---------------------------------------

columns(org.Rn.eg.db)

# Converting the gene IDs to gene names 
convertIDs <- function(ids, from, to, db, ifMultiple = c("putNA", "useFirst")) {
  stopifnot(inherits(db, "AnnotationDb"))
  ifMultiple <- match.arg(ifMultiple)
  suppressWarnings(selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to)))
  if (ifMultiple == "putNA") {
    duplicatedIds <- selRes[duplicated(selRes[,1]), 1]
    selRes <- selRes[! selRes[,1] %in% duplicatedIds,]
  }
  return(selRes[match(ids, selRes[,1] ), 2])
}

res$GeneID <- row.names(res)
res$gene_symbol <- convertIDs(row.names(res), "ENSEMBL", "SYMBOL", org.Rn.eg.db)

summary(res)

res_df <- as.data.frame(res)


#---------------------------------------
# Differentially expressed genes
#---------------------------------------

res_df$regulation_level <- "UNCHANGED"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
res_df$regulation_level[res_df$log2FoldChange > 0.5 & res_df$padj < 0.05] <- "UP"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
res_df$regulation_level[res_df$log2FoldChange < -0.5 & res_df$padj < 0.05] <- "DOWN"

write.table(res_df,
            file = "DESeq2_res.csv",
            sep = ",",
            row.names = F,
            col.names = T,
            quote = F)

res_df <- read.csv("DESeq2_res.csv", header = T)

res_df$regulation_level <- factor(res_df$regulation_level, levels = c("UP", "DOWN", "UNCHANGED"))

res_df <- res_df[!is.na(res_df$padj),]

pdf("volcano_plot.pdf")
ggplot(data=res_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = regulation_level),
       col = regulation_level,
       label = delabel) +
  geom_point(size = 1) +
  scale_color_manual(values = c("darkorchid", "royalblue", "lightgrey")) +
  xlab("Log2-Fold Change") +
  ylab("-log10(adjusted p-value)") +
  labs(color = "Regulation level") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# More visualizations
res_df_sorted <- res_df[order(res_df$padj),]

top.20 <- assay(rld)[row.names(assay(rld)) %in% res_df_sorted$GeneID[1:20],]
row.names(top.20) <- res_df_sorted[1:20,"gene_symbol"]

top.20.reshaped <- melt(top.20)
top.20.reshaped$Var3 <- rep(c("Neuropil", "Somata"), each = 60)

g1 <- ggplot(data = top.20.reshaped, aes(x = Var1, y = value, group = Var2)) +
  geom_point(aes(x = Var1, y = value, color = Var2), size = 2) +
  xlab("Genes") +
  ylab("Normalized read counts") +
  labs(color = "Samples") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

g2 <- ggplot(data = top.20.reshaped, aes(x = Var1, y = value, group = Var3)) +
  geom_point(aes(x = Var1, y = value, color = Var3), size = 2) +
  xlab("Genes") +
  ylab("Normalized read counts") +
  labs(color = "Samples") +
  ggtitle("Top 20 Significant DE Genes") +
  scale_color_manual(values = c("darkorchid", "royalblue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(g1, g2, ncol=1, nrow=2)