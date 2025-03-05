
# Load required packages
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

# Set environment
invisible(utils::memory.limit(100000))
library("BiocParallel")
register(SnowParam(6))
Sys.setenv("VROOM_CONNECTION_SIZE"=88888888)

# Set working directory
setwd("/home/ssloveff/lucas")

# Read data from cellranger output
Young_TN <- Read10X(data.dir = "Data/Young-TN/")
Young_3D <- Read10X(data.dir = "Data/Young-3d/")
Young_14D <- Read10X(data.dir = "Data/Young-14d/")
Aged_TN <- Read10X(data.dir = "Data/Aged-TN/")
Aged_3D <- Read10X(data.dir = "Data/Aged-3d/")
Aged_14D <- Read10X(data.dir = "Data/Aged-14d/")

# Create Seurat objects
Young_TN <- CreateSeuratObject(counts = Young_TN, project = "Young_TN", min.cells = 3, min.features = 200)
Young_3D <- CreateSeuratObject(counts = Young_3D, project = "Young_3D", min.cells = 3, min.features = 200)
Young_14D <- CreateSeuratObject(counts = Young_14D, project = "Young_14D", min.cells = 3, min.features = 200)
Aged_TN <- CreateSeuratObject(counts = Aged_TN, project = "Aged_TN", min.cells = 3, min.features = 200)
Aged_3D <- CreateSeuratObject(counts = Aged_3D, project = "Aged_3D", min.cells = 3, min.features = 200)
Aged_14D <- CreateSeuratObject(counts = Aged_14D, project = "Aged_14D", min.cells = 3, min.features = 200)

# Calculate mitochondrial gene percentage
Young_TN[["percent.mt"]] <- PercentageFeatureSet(Young_TN, pattern = "^mt-")
Young_3D[["percent.mt"]] <- PercentageFeatureSet(Young_3D, pattern = "^mt-")
Young_14D[["percent.mt"]] <- PercentageFeatureSet(Young_14D, pattern = "^mt-")
Aged_TN[["percent.mt"]] <- PercentageFeatureSet(Aged_TN, pattern = "^mt-")
Aged_3D[["percent.mt"]] <- PercentageFeatureSet(Aged_3D, pattern = "^mt-")
Aged_14D[["percent.mt"]] <- PercentageFeatureSet(Aged_14D, pattern = "^mt-")

# Quality control filtering
Young_TN <- subset(Young_TN, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
Young_3D <- subset(Young_3D, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
Young_14D <- subset(Young_14D, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
Aged_TN <- subset(Aged_TN, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
Aged_3D <- subset(Aged_3D, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
Aged_14D <- subset(Aged_14D, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

# Function for normalization and finding variable features
myfunction1 <- function(xxxx){
  xxxx <- NormalizeData(xxxx, normalization.method = "LogNormalize", scale.factor = 10000)
  xxxx <- FindVariableFeatures(xxxx, selection.method = "vst", nfeatures = 2000)
  return(xxxx)   
}

# Apply normalization to each sample
Young_TN <- myfunction1(Young_TN)
Young_3D <- myfunction1(Young_3D)
Young_14D <- myfunction1(Young_14D)
Aged_TN <- myfunction1(Aged_TN)
Aged_3D <- myfunction1(Aged_3D)
Aged_14D <- myfunction1(Aged_14D)

# Data integration
AB.anchors <- FindIntegrationAnchors(object.list = list(Young_TN,Young_3D,Young_14D,Aged_TN,Aged_3D,Aged_14D), 
                                    anchor.features = 2000,
                                    dims = 1:30)
AB <- IntegrateData(anchorset = AB.anchors, dims = 1:30)

# Data processing and dimensionality reduction
DefaultAssay(AB) <- "integrated"
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
AB <- FindNeighbors(AB, dims = 1:20)
AB <- FindClusters(AB, resolution = 0.4)
AB <- RunUMAP(AB, dims = 1:20)

# Add group information
AB$Age <- str_replace(AB$orig.ident,"_.*","")
AB$Group <- str_replace(AB$orig.ident,".*_","")
AB$Sample <- AB$orig.ident

# Set factor levels
AB@meta.data[["Sample"]] <- factor(AB@meta.data[["Sample"]], 
                                  levels=c("Young_TN","Young_3D","Young_14D","Aged_TN","Aged_3D","Aged_14D"))
AB@meta.data[["Group"]] <- factor(AB@meta.data[["Group"]], levels=c("TN","3D","14D"))
AB@meta.data[["Age"]] <- factor(AB@meta.data[["Age"]], levels=c("Young","Aged"))

# Save processed data
DefaultAssay(AB) <- "RNA"
save(AB, file="AB_new.RData")

# Plot UMAP visualization
p7 <- DimPlot(AB, reduction = "umap", group.by = "Age", pt.size=0.01,raster=FALSE) +
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "ident", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE) +
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2_Cluster.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')

# Find marker genes
DefaultAssay(AB) <- "RNA"
AB100 <- subset(AB, downsample=100)
markers100 <- FindAllMarkers(AB100, logfc.threshold = 0.25, min.pct = 0.5, only.pos = F, test.use = "wilcox")
write.table(markers100, file="markers100.txt", quote=F, sep="\t", row.names=F, col.names=T)

jjplot <- jjDotPlot(object = AB,
                    gene = c('Dpp4','Wnt2','Bmp7','Pi16','Gpr1','Icam1','Pparg','Fabp4','Cd36','F3','Fmo2','Clec11a','Wnt6','Sfrp5','Adipoq','Plin1','Car3','Spp1','Pecam1','Pdgfrb','Cspg4','Mcam','Mpz','Ptprc',"Hba-a1","Hba-a2","Plp1","Mbp",'Cd44','Itga6','Itgb1','Lgals3','Rhoa','Sp3','Bax','Fmod','Fli1','Myh11','Ebf2','Pdgfra','Cnn1','Myl9','Myocd','Tagln','Acta2','Fgf1','Ramp1','Fn1','Cd55','Dlk1','Col18a1','Mif','Flt1','Fth1','Postn','Tmsb4x','Dmkn','Krtdap','Apoe','Cd81','Ntrk3','Sparcl1','Thbs4','Tnc','Maoa','Rara','Acaca','Acly','Pnpla3','Lpl','Adam12','Bmper','Prdm16','Atxn1','Sfrp2','Igf1','Npr3','Nnat','Lep','Fgf14','Ces1f','Gsta3','Adrb3','Igfbp3','Cdkn1a','Penk','Pcsk1','Grn','Sort1','Gas6','Axl','Tspan5','Tmem219','Wnt5a','Lrpap1','Ackr2','Ryk','Egfr','Ngf','Cxcl12','Cst3','Timp1','Foxj3','Hsf1','Rgs5','Dcn','Lum','Col1a1','Apod','Mrtfa','Itga5','Apoc1','Cd9','Gas1','Epha3','Fabp5',"Cxcl14",'Gdf10',"Tmem219","Lrp1"),
                    gene.order = c('Dpp4','Wnt2','Bmp7','Pi16','Gpr1','Icam1','Pparg','Fabp4','Cd36','F3','Fmo2','Clec11a','Wnt6','Sfrp5','Adipoq','Plin1','Car3','Spp1','Pecam1','Pdgfrb','Cspg4','Mcam','Mpz','Ptprc',"Hba-a1","Hba-a2","Plp1","Mbp",'Cd44','Itga6','Itgb1','Lgals3','Rhoa','Sp3','Bax','Fmod','Fli1','Myh11','Ebf2','Pdgfra','Cnn1','Myl9','Myocd','Tagln','Acta2','Fgf1','Ramp1','Fn1','Cd55','Dlk1','Col18a1','Mif','Flt1','Fth1','Postn','Tmsb4x','Dmkn','Krtdap','Apoe','Cd81','Ntrk3','Sparcl1','Thbs4','Tnc','Maoa','Rara','Acaca','Acly','Pnpla3','Lpl','Adam12','Bmper','Prdm16','Atxn1','Sfrp2','Igf1','Npr3','Nnat','Lep','Fgf14','Ces1f','Gsta3','Adrb3','Igfbp3','Cdkn1a','Penk','Pcsk1','Grn','Sort1','Gas6','Axl','Tspan5','Tmem219','Wnt5a',"Lrpap1",'Ackr2','Ryk','Egfr','Ngf','Cxcl12','Cst3','Timp1','Foxj3','Hsf1','Rgs5','Dcn','Lum','Col1a1','Apod','Mrtfa','Itga5','Apoc1','Cd9','Gas1','Epha3','Fabp5',"Cxcl14",'Gdf10',"Tmem219","Lrp1"),
                    id = 'celltype',
                    cluster.order = c("Progenitor",		"Aregs",	"Preadipocyte",		"SMC",	"APC",	"Spp1+cells",	"Endo"),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)#+ggplot2:::coord_flip()
jjplot

ggsave(filename = "jjDotPlot celltype .pdf", plot = jjplot, device = 'pdf', width = 25, height = 11, units = 'cm')
