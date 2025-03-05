library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
library(ggpubr)
library(ggalluvial)
library(tidyverse)
invisible(utils::memory.limit(100000))
library("BiocParallel")
register(SnowParam(6))
Sys.setenv("VROOM_CONNECTION_SIZE"=88888888)

setwd("G:/Lucas/cellchat")

Aregs <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("Aregs")]
SMC <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("SMC")]
Endo <- AB[,AB@meta.data[["celltype"]] %in% c("Endo")]
AB <- merge(Aregs,SMC)

load("AB.RData")
rm("AB")
#
AB$age=str_replace(AB$age,"Sham3d","Sham")
AB$age=str_replace(AB$age,"Sham7d","Sham")

AB@meta.data[["age"]]<-factor(AB@meta.data[["age"]], levels=c('Young_TN','Young_3D','Young_14D','Aged_TN','Aged_3D','Aged_14D'))

age.list<-SplitObject(AB, split.by = "Sample")
#
data.input  <- age.list[["Aged_14D"]]@assays[["RNA"]]@data
celltype  <- age.list[["Aged_14D"]]@meta.data[["celltype"]]
data.input[1:4,1:4]
identity = data.frame(group = age.list[["Aged_14D"]]$celltype, row.names = names(age.list[["Aged_14D"]]$celltype)) 
head(identity)
unique(identity$group) 
table(identity$group)

meta <- data.frame(labels = age.list[["Aged_14D"]]$celltype, row.names = names(identify))
cellchat <- createCellChat(object = data.input,meta = meta, group.by = "labels")
cellchat
summary(cellchat)

cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
head(cellchat@meta)
## set "labels" as default cell identity
cellchat <- setIdent(cellchat, ident.use = "labels") 

Young_TN <- cellchat
Young_3D <- cellchat
Young_14D <- cellchat
Aged_TN <- cellchat
Aged_3D <- cellchat
Aged_14D <- cellchat

cellchat <- Young_TN
cellchat <- Young_3D
cellchat <- Young_14D
cellchat <- Aged_TN
cellchat <- Aged_3D
cellchat <- Aged_14D

CellChatDB <- CellChatDB.mouse 
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
# set the used database in the object
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") #"Secreted Signaling" "ECM-Receptor" "Cell-Cell Contact"
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 
######plan("multiprocess", workers = 4) # do parallel  
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.mouse)
##

cellchat <- computeCommunProb(cellchat,type = "triMean", population.size = F)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

Young_TN <- cellchat
Young_3D <- cellchat
Young_14D <- cellchat
Aged_TN <- cellchat
Aged_3D <- cellchat
Aged_14D <- cellchat

saveRDS(Young_TN,"Young_TN.rds")
saveRDS(Young_3D,"Young_3D.rds")
saveRDS(Young_14D,"Young_14D.rds")
saveRDS(Aged_TN,"Aged_TN.rds")
saveRDS(Aged_3D,"Aged_3D.rds")
saveRDS(Aged_14D,"Aged_14D.rds")

load("Sham.rds")
load("Young_3D.rds")
load("Young_14D.rds")
load("Aged_TN.rds")
load("Aged_3D.rds")
load("Aged_14D.rds")

con.list <- list(Young_TN=Young_TN,Young_3D=Young_3D,Young_14D=Young_14D,Aged_TN=Aged_TN,Aged_3D=Aged_3D,Aged_14D=Aged_14D)

cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)

con.list <- list(Young_TN=Young_TN,Young_3D=Young_3D,Young_14D=Young_14D,Aged_TN=Aged_TN,Aged_3D=Aged_3D,Aged_14D=Aged_14D
)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)

levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(1,2), comparison = c(1,2,3,4,5,6), angle.x = 45)
p
ggsave("Aregs-SMC(ligand)Compare_LR.pdf", p, width = 5, height = 13)
