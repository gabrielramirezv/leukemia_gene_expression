library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

genes <- read.table("lymphoma_Burkitt_sobreexpresados_high_vs_low.txt")
genes <- genes$V1


annotations <- AnnotationDbi::select(org.Hs.eg.db, 
                                     keys = genes, 
                                     columns = c("ENSEMBL", 
                                                 "SYMBOL", 
                                                 "ENTREZID", 
                                                 "GO"), 
                                     keytype = "ENSEMBL")

ego <- enrichGO(gene = genes,
                OrgDb = "org.Hs.eg.db",
                keyType = "ENSEMBL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

barplot(ego, showCategory = 10) + 
  ggtitle("Enrichment Analysis of lymphoma Burkitt Genes") +
  theme_classic()

gseaKEGG <- gseKEGG(geneList = genes,
                organism = "hsa",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                minGSSize = 10)
