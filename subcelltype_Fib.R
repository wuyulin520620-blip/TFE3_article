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
immune.combined = subset(sce,Celltype %in% c('Fibroblast'))
DefaultAssay(immune.combined) = 'RNA'
sceList <- SplitObject(immune.combined, split.by = "orig.ident")
single.ob=merge(sceList[[1]],sceList[2:length(sceList)])
DefaultAssay(single.ob) = 'RNA'
scRNA_harmony = single.ob
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
##整合
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:20)
testAB.integrated <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution =  0.1)

Idents(testAB.integrated) = 'seurat_clusters'
DefaultAssay(testAB.integrated) = 'RNA'
sce.markers <- FindAllMarkers(object = testAB.integrated, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
library(dplyr) 
# 不同seurat版本的 avg_logFC 不一样 
top5 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

SCE = testAB.integrated

SCE = RenameIdents(SCE,
                   '0' = 'Collagen_Fib',
                   '1' = 'Collagen_Fib',
                   '4' = 'Collagen_Fib',
                   '5' = 'Collagen_Fib',
                   '6' = 'Collagen_Fib',
                   '7' = 'NACM1_Fib')

SCE <- AddMetaData(SCE,SCE@active.ident,col.name = "Celltype")
saveRDS(SCE,'./scRNA_Fid.rds')


