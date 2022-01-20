#---------------------------------------
# Gene Ontology Analysis
#
#
# This script comprises hands-on code
#---------------------------------------

library(topGO)
library(ggplot2)

# The script follows the "differential expression analysis" part
# Data preparation
res_df_sorted <- res_df[order(res_df$padj),]

genes_name <- row.names(res_df_sorted)

upRegulated.genes <- which(res_df$padj < 0.05 & res_df$log2FoldChange > 0)
downRegulated.genes <- which(res_df$padj < 0.05 & res_df$log2FoldChange < 0)

all_genes_names <- rownames(res_df)

genes_up <- rownames(res_df)[upRegulated.genes]
genes_down <- rownames(res_df)[downRegulated.genes]

genelist_up <- factor(as.integer(all_genes_names %in% genes_up))
names(genelist_up) <- all_genes_names

genelist_down <- factor(as.integer(all_genes_names %in% genes_down))
names(genelist_down) <- all_genes_names

# Provided function which maps gene identifiers to GO terms
allGO2genes <- annFUN.org(whichOnto = "ALL",
                          feasibleGenes = NULL,
                          mapping = "org.Rn.eg.db",
                          ID = "ensembl")

getGO.infos <- function(ontology_type, geneList_type) {
  GO.object <- new("topGOdata",
                   ontology = ontology_type,
                   allGenes = geneList_type,
                   annot = annFUN.GO2genes,
                   GO2genes = allGO2genes,
                   nodeSize = 10)
  return(GO.object)
}

GO.Up.BP <- getGO.infos("BP", genelist_up)
GO.Up.MF <- getGO.infos("MF", genelist_up)
GO.Up.CC <- getGO.infos("CC", genelist_up)

GO.Down.BP <- getGO.infos("BP", genelist_down)
GO.Down.MF <- getGO.infos("MF", genelist_down)
GO.Down.CC <- getGO.infos("CC", genelist_down)

# Enrichment tests using Fisher's exact test
resultFisher <- function(x) {
  test <- runTest(x, statistic = "fisher")
  return(test)
}

result.Up.BP <- resultFisher(GO.Up.BP)
result.Up.MF <- resultFisher(GO.Up.MF)
result.Up.CC <- resultFisher(GO.Up.CC)

result.Down.BP <- resultFisher(GO.Down.BP)
result.Down.MF <- resultFisher(GO.Down.MF)
result.Down.CC <- resultFisher(GO.Down.CC)

# Parsing the table
parse_tables <- function(GO_data, statistics)
{
  goEnrichment <- GenTable(GO_data, weightFisher = statistics, topNodes = 20)
  sub("< ", "", goEnrichment$weightFisher)
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels = rev(goEnrichment$Term))
  goEnrichment$weightFisher <- as.numeric(sub("< ", "", goEnrichment$weightFisher))  
  goEnrichment
}

GOres_up_bp <- parse_tables(GO.Up.BP, result.Up.BP)
GOres_up_mf <- parse_tables(GO.Up.MF, result.Up.MF)
GOres_up_cc <- parse_tables(GO.Up.CC, result.Up.CC)

GOres_down_bp <- parse_tables(GO.Down.BP, result.Down.BP)
GOres_down_mf <- parse_tables(GO.Down.MF, result.Down.MF)
GOres_down_cc <- parse_tables(GO.Down.CC, result.Down.CC)

plot_GO <- function(GO_data, Ontology, Regulation, use_color) {
  GO_data$log_weightFisher <- (- log10(as.numeric(GO_data$weightFisher)))
  ggplot(GO_data, 
         aes(x = GO_data$log_weightFisher,
             y = GO_data$Term)) +
    geom_point(aes(size = GO_data$Significant),
               colour = use_color) +
    scale_size_area(name = "Gene counts") +
    xlab("Enrichment (- log10 Pvalue)") +
    ylab(Ontology) +
    ggtitle(Regulation) +
    scale_x_continuous() +
    theme_bw() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black", size = 7),
      axis.title=element_text(size=12),)
}

p1 <- plot_GO(GOres_up_bp, "Biological Proccess", "TopGO Up (fisher's exact test)", "darkorchid")
p2 <- plot_GO(GOres_up_mf, "Molecular Function", "TopGO Up (fisher's exact test)", "darkorchid")
p3 <- plot_GO(GOres_up_cc, "Cellular Component", "TopGO Up (fisher's exact test)", "darkorchid")

pdf("All_types_upregulated.pdf")
grid.arrange(grobs = list(p1,p2,p3))
dev.off()

p4 <- plot_GO(GOres_down_bp, "Biological Proccess", "TopGO Down (fisher's exact test)", "royalblue")
p5 <- plot_GO(GOres_down_mf, "Molecular Function", "TopGO Down (fisher's exact test)", "royalblue")
p6 <- plot_GO(GOres_down_cc, "Cellular Component", "TopGO Down (fisher's exact test)", "royalblue")

pdf("All_types_downregulated.pdf")
grid.arrange(grobs = list(p4,p5,p6))
dev.off()