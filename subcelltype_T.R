library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(scRNAtoolVis)
library(scales)
library(harmony)

scRNA = readRDS('./scRNA.rds')
immune.combined = subset(scRNA,Celltype %in% c('T cell'))

immune.combined <- FindVariableFeatures(object = immune.combined,selection.method = 'vst', nfeatures = 2000)
ifnb.list <- SplitObject(immune.combined, split.by = "orig.ident")
testAB.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:50)
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:50)
DefaultAssay(testAB.integrated) <- "integrated"
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:50)
testAB.integrated <- RunTSNE(testAB.integrated, dims = 1:50)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:50)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.4)


DefaultAssay(SCE) = 'RNA'
Idents(SCE) = 'seurat_clusters'
DefaultAssay(SCE) = 'RNA'
sce.markers <- FindAllMarkers(object = SCE, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
library(dplyr) 
# 不同seurat版本的 avg_logFC 不一样 
top5 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
testAB.integrated <- CellCycleScoring(testAB.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

SCE = testAB.integrated

celltype = data.frame(
  cluster = c(0,1,3,4,5,6,8),
  celltype = c('Cd8_Tex_TOX','Memory T cells','Tn',
               'Treg','NK','NK','Cd8_Tex_TOX'
  )
)
for(i in 1:nrow(celltype)){
  SCE@meta.data[which(SCE@meta.data$seurat_clusters == celltype$cluster[i]),'Celltype'] <- celltype$celltype[i]}

saveRDS(SCE,'./scRNA_T.rds')



