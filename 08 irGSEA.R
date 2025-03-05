
library(Seurat)
library(ggplot2)

library(SeuratData)


scRNA <- scRNA
scRNA <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("Aregs")]
Idents(scRNA) <- scRNA$Sample

# plot
DimPlot(scRNA, reduction = "umap",
        group.by = "celltype",label = T) + NoLegend()
library(AUCell)
library(UCell)
library(irGSEA)
library(xlsx)
library(readxl)

Idents(scRNA) <- scRNA$celltype

msigdbr::msigdbr_species()
msigdbr::msigdbr_collections()
###
markers.to.plot <- c("Calm1",	"Calm2",	"Calm3",	"Fkbp1a",	"Nfatc1",	"Nfatc2",	"Nfatc3",	"Ppp3ca",	"Ppp3cb",	"Ppp3r1")

dp2 <- DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2, dot.scale = 9) + RotatedAxis() #+ggplot2:::coord_flip()#####seurat_clusters

dp2
####
geneset <- read_xlsx("Gene.xlsx",col_names = T)


scRNA <- irGSEA.score(object = scRNA, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = T, geneset = geneset, msigdb = T, 
                             species = "Rattus norvegicus", category = "GO",
                             subcategory = 'GO:BP', geneid = "symbol",
                             method = c("AUCell"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

Seurat::Assays(scRNA)

result.dge <- irGSEA.integrate(object = scRNA, 
                               group.by = "celltype",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell", "UCell", "singscore", 
                                          "ssgsea","JASMINE"))
class(result.dge)

densityheatmap <- irGSEA.densityheatmap(object = scRNA,
                                        method = "AUCell",
                                        show.geneset = "Gene")
densityheatmap
ggsave(filename = "densityheatmap-TGFB+P53.pdf", plot = densityheatmap, device = 'pdf', width = 15, height = 15, units = 'cm')











