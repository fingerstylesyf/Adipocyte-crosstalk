library(Matrix)
library(Seurat)
library(ROGUE)

AB <- AB[,AB@meta.data[["celltype"]] %in% c("Aregs")]

AB$UMAP_1 <- AB@reductions$umap@cell.embeddings[,1]
AB$UMAP_2 <- AB@reductions$umap@cell.embeddings[,2]

AB$barcode <- tmp1
write.csv(AB@meta.data, 
          file='metadata.csv', quote=F, row.names=F)

library(Matrix)
counts_matrix <- GetAssayData(AB, assay='RNA', slot='counts')
expr <- data.frame(counts_matrix)
write.csv(p1, file = "counts_matrix.csv", row.names = TRUE)

expr[1:5, 1:4]

meta <- read.csv("metadata.csv")
head(meta)

expr <- matr.filter(expr, min.cells = 10, min.genes = 10)

ent.res <- SE_fun(expr)
head(ent.res)

SEplot(ent.res)

rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value

rogue.res <- rogue(expr, labels = meta$celltype, samples = meta$Sample, platform = "UMI", span = 0.6)
rogue.res

rogue.boxplot(rogue.res)
write.csv(rogue.res, file = "WE.rogue.res.csv", row.names = TRUE)












