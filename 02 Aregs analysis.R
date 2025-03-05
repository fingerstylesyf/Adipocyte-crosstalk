library(Matrix)
library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(DoubletFinder)
library(magrittr)
library(stringr)
library(scRNAtoolVis)

invisible(utils::memory.limit(100000))
library("BiocParallel")
register(SnowParam(6))
Sys.setenv("VROOM_CONNECTION_SIZE"=88888888)

setwd("/home/ssloveff/lucas/02_subtype/Aregs")
getwd()

# Data normalization and feature selection
normalize_and_find_features <- function(data) {
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  return(data)
}

scRNA <- normalize_and_find_features(scRNA)

# Data integration
object.list <- SplitObject(scRNA, split.by = "orig.ident")
scRNA.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = 2000, dims = 1:30)
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)

# Standard workflow for visualization and clustering
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
ElbowPlot(scRNA, ndims = 30)
save(scRNA, file="scRNA.RData")

# Clustering
scRNA <- FindNeighbors(scRNA, dims = 1:12)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims = 1:12)
DimPlot(scRNA, reduction = "umap")

# Extract cell counts
table(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$seurat_clusters)

# Reordering groups
scRNA@meta.data[["Sample"]] <- factor(scRNA@meta.data[["Sample"]], levels=c("Young_TN","Young_3D","Young_14D","Aged_TN","Aged_3D","Aged_14D"))
scRNA@meta.data[["Group"]] <- factor(scRNA@meta.data[["Group"]], levels=c("TN","3D","14D"))
scRNA@meta.data[["Age"]] <- factor(scRNA@meta.data[["Age"]], levels=c("Young","Aged"))

DefaultAssay(scRNA) <- "RNA"
save(scRNA, file="scRNA_new.RData")

# UMAP visualization
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Age", pt.size=0.01, raster=FALSE) +
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p2 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", pt.size=0.01, label = TRUE, repel = TRUE, raster=FALSE) +
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
plot_grid(p1, p2, align = "v", ncol = 2)
ggsave(filename = "umap2_Cluster.pdf", plot = last_plot(), device = 'pdf', width = 43, height = 17, units = 'cm')

# Finding marker genes
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = FALSE, test.use = "wilcox")
write.table(markers, file="markers.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Removing unwanted cells
scRNA <- scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("0","1","2","3","4","5","6","7")]

# Second clustering and dimensionality reduction
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
ElbowPlot(scRNA, ndims = 30)
save(scRNA, file="scRNA3.RData")

scRNA <- FindNeighbors(scRNA, dims = 1:10)
scRNA <- FindClusters(scRNA, resolution = 0.2)
scRNA <- RunUMAP(scRNA, dims = 1:10)
DimPlot(scRNA, reduction = "umap")

# UMAP visualization
p3 <- DimPlot(scRNA, reduction = "umap", group.by = "age", pt.size=0.1)
p4 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", pt.size=0.01, label = TRUE, repel = TRUE)
plot_grid(p3, p4, align = "v", ncol = 2)
ggsave(filename = "umap2.pdf", plot = last_plot(), device = 'pdf', width = 34, height = 14.5, units = 'cm')

# Finding marker genes on a subset
scRNA100 <- subset(scRNA, downsample=400)
markers <- FindAllMarkers(scRNA100, logfc.threshold = 0.25, min.pct = 0.25, only.pos = FALSE, test.use = "wilcox")
write.table(markers, file="markers-cluster.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Assigning new cell type identities
current.cluster.ids <- c("0","1","2","3")
new.cluster.ids <- c("Preadipocyte", "Progenitor", "Aregs", "Progenitor")
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$celltype)

Idents(scRNA) <- factor(scRNA$celltype, levels = c("Preadipocyte", "Progenitor", "Aregs"))
DimPlot(scRNA, reduction = "umap", label = TRUE)
save(scRNA, file="Aregs_renamed.RData")

# Finding marker genes again
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = FALSE, test.use = "wilcox")
write.table(markers, file="markers-celltype.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

