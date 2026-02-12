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
immune.combined = subset(scRNA,Celltype %in% c('Endothelial cell'))

DefaultAssay(immune.combined) = 'RNA'
sceList <- SplitObject(immune.combined, split.by = "orig.ident")
single.ob=merge(sceList[[1]],sceList[2:length(sceList)])
DefaultAssay(single.ob) = 'RNA'
scRNA_harmony = single.ob
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
##整合
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:30)
testAB.integrated <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution =  0.12)

DefaultAssay(testAB.integrated) = 'RNA'
Idents(testAB.integrated) = 'seurat_clusters'
sce.markers <- FindAllMarkers(object = testAB.integrated, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
# 不同seurat版本的 avg_logFC 不一样 
top5 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
testAB.integrated <- CellCycleScoring(testAB.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

celltype = data.frame(
  cluster = c(0,2,3,4),
  celltype = c('VWF_EC','NTN4_EC','IGFBP3_EC',
               'Lymp_EC'
  )
)
for(i in 1:nrow(celltype)){
  SCE@meta.data[which(SCE@meta.data$seurat_clusters == celltype$cluster[i]),'Celltype'] <- celltype$celltype[i]}

saveRDS(SCE,'./scRNA_Endo.rds')


