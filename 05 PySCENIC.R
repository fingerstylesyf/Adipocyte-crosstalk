
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)

#
sce_SCENIC <- open_loom("sce_SCENIC.loom")

regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')

cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"
sub_regulonAUC <- regulonAUC

rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])

rssPlot <- plotRSS(
    rss,
    order_rows = TRUE,
    thr=0.1,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')
rssPlot

rss_data <- rssPlot$plot$data

rss_data<-dcast(rss_data, 
                Topic~rss_data$cellType,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
colnames(rss_data)
col_ann <- data.frame(group= c(rep("SMC",1),
                               rep("Endo",1),
                               rep("Aregs",1),
                               rep("Progenitor",1),
                               rep("APC",1),
                               rep("Spp1+cells",1),
                               rep("Preadipocyte",1)))#
rownames(col_ann) <- colnames(rss_data)

col <- list(group=groupcol)

text_columns <- sample(colnames(rss_data),0)

write.csv(rss_data, 'Celltype-TFs.csv', quote=F) 
