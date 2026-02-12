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
immune.combined = subset(scRNA,Celltype %in% c('Neutrophil','Macrophage','Mast','Prof'))
immune.combined <- FindVariableFeatures(object = immune.combined,selection.method = 'vst', nfeatures = 3000)
ifnb.list <- SplitObject(immune.combined, split.by = "orig.ident")
testAB.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:30)
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:30)
DefaultAssay(testAB.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 30, verbose = FALSE)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:30)
testAB.integrated <- RunTSNE(testAB.integrated, dims = 1:30)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:30)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.3)

s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
testAB.integrated <- CellCycleScoring(testAB.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(testAB.integrated) = 'RNA'
Idents(testAB.integrated) = 'seurat_clusters'
sce.markers <- FindAllMarkers(object = testAB.integrated, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
library(dplyr) 
# 不同seurat版本的 avg_logFC 不一样 
top5 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
SCE = testAB.integrated 
SCE = RenameIdents(SCE,
                   '0' = 'MRC1_RTM',
                   '1' = 'GPNMB_TAM',
                   '2' = 'CLEC4F_RTM',
                   '3' = 'Epi',
                   '4' = 'Prof.',
                   '5' = 'GPNMB_TAM',
                   '6' = 'Inter.Monocyte',
                   '7' = 'Fib',
                   '8' = 'Neutrophils',
                   '9' = 'Mast',
                   '10' = 'DC',
                   '11' = 'DC',
                   '12' = 'MRC1_RTM',
                   '13' = 'DC')
SCE <- AddMetaData(SCE,SCE@active.ident,col.name = "Celltype_Mac")
saveRDS(SCE,'./scRNA_Mac.rds')










