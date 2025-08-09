library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(pathview)

dea.res <- read.table("DESeq2_results_AML_high_vs_low.txt", header = TRUE, sep = "\t", row.names = 6)
dea.res.genes <- rownames(dea.res)

annotations <- AnnotationDbi::select(org.Hs.eg.db, 
                                     keys = dea.res.genes, 
                                     columns = c("ENSEMBL", 
                                                 "SYMBOL", 
                                                 "ENTREZID", 
                                                 "GO"), 
                                     keytype = "ENSEMBL")

annotations <- annotations[!duplicated(annotations$ENSEMBL), ]

rownames(annotations) <- annotations$ENSEMBL

dea.res <- merge(dea.res, annotations, by = "row.names", all = TRUE)

rownames(dea.res) <- dea.res$ENSEMBL

entrez.ids <- dea.res$ENTREZID
names(entrez.ids) <- rownames(dea.res)
entrez.ids <- entrez.ids[!duplicated(entrez.ids)]

dea.res.lfc <- dea.res$log2FoldChange
names(dea.res.lfc) <- dea.res$ENTREZID
dea.res.lfc <- dea.res.lfc[entrez.ids]

dea.res.lfc <- sort(dea.res.lfc, decreasing = TRUE)

gseaKEGG <- gseKEGG(geneList = dea.res.lfc,
                organism = "hsa",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                minGSSize = 10)
gseaKEGG.res <- gseaKEGG@result

gseaplot(gseaKEGG, geneSetID = "hsa04382", title = "Pathways in cancer")
pathview(gene.data = dea.res.lfc, pathway.id = "hsa04382", species = "hsa",
         limit = list(gene = 2, cpd = 1))
