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
seuratdata = subset(scRNA,Celltype == 'Epithelial')

sceList <- SplitObject(seuratdata, split.by = "orig.ident")
single.ob=merge(sceList[[1]],sceList[2:length(sceList)])
DefaultAssay(single.ob) = 'RNA'
scRNA_harmony <- NormalizeData(single.ob) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)

##整合
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA.integrated <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution =  0.3)

DefaultAssay(scRNA.integrated) = 'RNA'
Idents(scRNA.integrated) = 'seurat_clusters'
DefaultAssay(scRNA.integrated) = 'RNA'
sce.markers <- FindAllMarkers(object = scRNA.integrated, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
library(dplyr) 
# 不同seurat版本的 avg_logFC 不一样 
top5 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

scRNA.integrated = RenameIdents(scRNA.integrated,
                   '0' = 'CDH11_epi',
                   '1' = 'SPP1_epi',
                   '2' = 'PTGDS_epi',
                   '3' = 'ITGA3_epi',
                   '4' = 'PT',
                   '5' = 'TAL',
                   '6' = 'PC',
                   '8' = 'TAL',
                   '9' = 'Type_B_Intercalated_cell',
                   '10' = 'CDH11_epi',
                   '11' = 'ITGA3_epi')
scRNA.integrated <- AddMetaData(scRNA.integrated,scRNA.integrated@active.ident,col.name = "Celltype")
saveRDS(scRNA.integrated,'./scRNA_Epi.rds')




