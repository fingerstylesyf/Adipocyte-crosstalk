# Loading required libraries
library(multtest)
library(Seurat)
library(dplyr)
library(mindr)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)
library(scRNAtoolVis)
library(Scillus)

# Set working directory
setwd("G:/Lucas/02 Subgroup Analysis/Adipo")

# Custom function for preprocessing
myfunction1 <- function(xxxx){
  xxxx <- NormalizeData(xxxx, normalization.method = "LogNormalize", scale.factor = 10000)  # Normalize data
  xxxx <- FindVariableFeatures(xxxx, selection.method = "vst", nfeatures = 2000)  # Find variable features (VST variance-stabilizing transformation)
  return(xxxx)   
}

# Apply custom function to the data
scRNA <- myfunction1(scRNA)

# Data integration
object.list <- SplitObject(scRNA, split.by = "orig.ident")  # Create new integration grouping
scRNA.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = 2000, dims = 1:30)  # Find integration anchors
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)  # Integrate data

# Set default assay and scale data
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))

# Perform PCA and visualize results
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
plot2 <- ElbowPlot(scRNA, ndims = 30)
plot2

# Save the data after PCA
save(scRNA, file = "scRNA.RData")

# Clustering and UMAP visualization
scRNA <- FindNeighbors(scRNA, dims = 1:11)
scRNA <- FindClusters(scRNA, resolution = 0.4)  # Adjust resolution for clustering
scRNA <- RunUMAP(scRNA, dims = 1:11)  # Perform UMAP

# Plot UMAP
p7 <- DimPlot(scRNA, reduction = "umap", group.by = "Age", pt.size = 0.1) +
  theme(plot.title = element_text(size = 0), strip.text = element_text(size = 20), axis.title = element_text(size = 20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5, label = TRUE, repel = TRUE) +
  theme(plot.title = element_text(size = 0), strip.text = element_text(size = 20), axis.title = element_text(size = 20))
umap1 <- plot_grid(p7, p8, align = "v", ncol = 2)
ggsave(filename = "umap.pdf", plot = umap1, device = 'pdf', width = 34, height = 14.5, units = 'cm')

# Extract cell counts and inspect data
table(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$seurat_clusters)
table(scRNA@meta.data$celltype)
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["Sample"]])

# Save data
DefaultAssay(scRNA) <- "RNA"
save(scRNA, file = "scRNA2.RData")

# Differential gene expression analysis
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = FALSE, test.use = "wilcox")
write.table(markers, file = "markers-cluster.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Dot plot for selected genes
jjplot <- jjDotPlot(object = scRNA,
                    gene = c('Dpp4', 'Wnt2', 'Bmp7', 'Pi16', 'Gpr1', 'Icam1', 'Pparg', 'Fabp4', 'Cd36', 'F3', 'Fmo2', 'Clec11a', 
                             'Wnt6', 'Sfrp5', 'Adipoq', 'Plin1', 'Car3', 'Spp1', 'Pecam1', 'Pdgfrb', 'Mpz', 'Ptprc', 'Hba-a1', 'Hba-a2',
                             'Plp1', 'Mbp', 'Cd44', 'Itga6', 'Itgb1', 'Lgals3', 'Rhoa', 'Sp3', 'Bax', 'Fmod', 'Fli1', 'Myh11', 'Ebf2', 
                             'Pdgfra', 'Cnn1', 'Myl9', 'Myocd', 'Tagln', 'Acta2', 'Fgf1', 'Ramp1', 'Fn1', 'Cd55', 'Dlk1', 'Col18a1', 
                             'Mif', 'Flt1', 'Fth1', 'Postn', 'Tmsb4x', 'Dmkn', 'Krtdap', 'Apoe', 'Cd81', 'Ntrk3', 'Sparcl1', 'Thbs4', 
                             'Tnc', 'Maoa', 'Rara', 'Acaca', 'Acly', 'Pnpla3', 'Lpl', 'Adam12', 'Bmper', 'Prdm16', 'Atxn1', 'Sfrp2', 'Igf1',
                             'Trpv1', 'Tnfrsf9', 'Id1', 'Ucp1', 'Gk', 'Esrrg'),
                    gene.order = c('Dpp4', 'Wnt2', 'Bmp7', 'Pi16', 'Gpr1', 'Icam1', 'Pparg', 'Fabp4', 'Cd36', 'F3', 'Fmo2', 'Clec11a', 
                                   'Wnt6', 'Sfrp5', 'Adipoq', 'Plin1', 'Car3', 'Spp1', 'Pecam1', 'Pdgfrb', 'Mpz', 'Ptprc', 'Hba-a1', 
                                   'Hba-a2', 'Plp1', 'Mbp', 'Cd44', 'Itga6', 'Itgb1', 'Lgals3', 'Rhoa', 'Sp3', 'Bax', 'Fmod', 'Fli1', 
                                   'Myh11', 'Ebf2', 'Pdgfra', 'Cnn1', 'Myl9', 'Myocd', 'Tagln', 'Acta2', 'Fgf1', 'Ramp1', 'Fn1', 'Cd55', 
                                   'Dlk1', 'Col18a1', 'Flt1', 'Fth1', 'Postn', 'Tmsb4x', 'Dmkn', 'Krtdap', 'Apoe', 'Cd81', 'Ntrk3', 
                                   'Sparcl1', 'Acaca', 'Acly', 'Pnpla3', 'Lpl', 'Adam12', 'Bmper', 'Prdm16', 'Ucp1', 'Gk', 'Esrrg'),
                    cluster.order = 0:6,
                    ytree = FALSE,
                    rescale = TRUE,
                    rescale.min = 0,
                    rescale.max = 2,
                    dot.max = 8)
jjplot
ggsave(filename = "jjDotPlot_after_clustering.pdf", plot = jjplot, device = 'pdf', width = 34, height = 13, units = 'cm')

# Remove unwanted clusters
scRNA <- scRNA[, scRNA@meta.data[["seurat_clusters"]] %in% c("0", "1", "2", "3", "4", "5")]

# Second round of PCA and UMAP
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
plot2 <- ElbowPlot(scRNA, ndims = 30)
plot2
save(scRNA, file = "scRNA3.RData")

# Clustering and UMAP visualization
scRNA <- FindNeighbors(scRNA, dims = 1:7)
scRNA <- FindClusters(scRNA, resolution = 0.2)
scRNA <- RunUMAP(scRNA, dims = 1:7)

# Plot UMAP again
DimPlot(scRNA, reduction = "umap")

# Save the processed data
DefaultAssay(scRNA) <- "RNA"
save(scRNA, file = "scRNA4.RData")

# Differential gene expression analysis
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = FALSE, test.use = "wilcox")
write.table(markers, file = "markers-cluster.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Reassign cell types based on clusters
current.cluster.ids <- c("0", "1", "2", "3", "4")
new.cluster.ids <- c("SMC", "SMC", "APC", "APC", "APC")
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)

# Check celltype distribution
table(scRNA@meta.data$celltype)

# Rename cell types
Idents(scRNA) <- factor(scRNA$celltype, levels = c("SMC", "APC"))

# Final UMAP visualization
DimPlot(scRNA, reduction = "umap", label = TRUE)

# Save the final data
save(scRNA, file = "SMC_celltype.RData")

# Differential gene expression analysis for renamed cell types
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = FALSE, test.use = "wilcox")
write.table(markers, file = "markers-celltype.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)