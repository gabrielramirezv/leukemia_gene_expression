# Import packages
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(pathview)

# Load the DESeq2 results
dea.res <- read.table("DESeq2_results_lymphoma_Hodgkin_high_vs_low.txt", 
                      header = TRUE, 
                      sep = "\t", 
                      row.names = 6)

# Ensure the row names are set correctly
dea.res.genes <- rownames(dea.res)

# Get annotations for the genes
annotations <- AnnotationDbi::select(org.Hs.eg.db, 
                                     keys = dea.res.genes, 
                                     columns = c("ENSEMBL", 
                                                 "SYMBOL", 
                                                 "ENTREZID", 
                                                 "GO"), 
                                     keytype = "ENSEMBL")

# Remove duplicates based on ENSEMBL IDs
annotations <- annotations[!duplicated(annotations$ENSEMBL), ]

# Ensure that the annotations are in the same order as the DEA results
rownames(annotations) <- annotations$ENSEMBL

# Merge DEA results with annotations
dea.res <- merge(dea.res, annotations, by = "row.names", all = TRUE)

# Set ENSEMBL IDs as row names
rownames(dea.res) <- dea.res$ENSEMBL

# Get ENTREZ IDs and remove duplicates
entrez.ids <- dea.res$ENTREZID
names(entrez.ids) <- rownames(dea.res)
entrez.ids <- entrez.ids[!duplicated(entrez.ids)]

# Get log2 fold changes for the ENTREZ IDs
dea.res.lfc <- dea.res$log2FoldChange
names(dea.res.lfc) <- dea.res$ENTREZID
dea.res.lfc <- dea.res.lfc[entrez.ids]

# Sort the log2 fold changes in decreasing order
dea.res.lfc <- sort(dea.res.lfc, decreasing = TRUE)

# Perform Gene Set Enrichment Analysis (GSEA) using KEGG pathways
gseaKEGG <- gseKEGG(geneList = dea.res.lfc,
                organism = "hsa",
                pvalueCutoff = 0.05,
                minGSSize = 20)

# Extract the results from the GSEA object
gseaKEGG@result

# Generate the GSEA plot for a specific pathway
gseaplot(gseaKEGG, geneSetID = "hsa05169", title = "Epstein-Barr virus infection")

# Visualize gene expression data for the specific pathway
pathview(gene.data = dea.res.lfc, pathway.id = "hsa05169", species = "hsa",
         limit = list(gene = 2, cpd = 1))
