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
library(GSVA)
library(org.Hs.eg.db)
library(GSEABase)
library(clusterProfiler)
library(DOSE)

scRNA = readRDS('./scRNA.rds')
l <- data.frame()
cells = unique(scRNA$Celltype)

for (i in cells){
  Idents(scRNA) = 'Celltype'
  cell = subset(scRNA,idents = i)
  DefaultAssay(cell) = 'RNA'
  Idents(cell) = 'Group'
  diffgene = FindMarkers(cell,ident.1 = 'Tumor',logfc.threshold = 0.25,
                         test.use = "wilcox",
                         min.pct = 0.15)
  diffgene$cluster = i
  diffgene$gene = rownames(diffgene)
  l <- rbind(l, diffgene)
}


expr <- AverageExpression(scRNA, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #过滤细胞表达量全为零的基因
expr <- as.matrix(expr)
ss = sort(as.vector(unique(scRNA$Celltype)) , decreasing = F)

human_KEGG_Set  <- getGmt('./h.all.v7.5.1.symbols.gmt') 
gsva.kegg <- gsva(expr, gset.idx.list = human_KEGG_Set,
                  kcdf="Gaussian",
                  method = "gsva",
                  parallel.sz=1)

a = c()
ds = as.data.frame(gsva.kegg) 
for (i in 1:dim(ds)[2] ){
  path = ds[ order (-ds[,i]),] %>%  head(50) %>%  rownames()
  a = unique(c(a,path))
}
gsva.kegg = gsva.kegg[a,]

rownames(gsva.kegg) = tolower( rownames(gsva.kegg) )
options(repr.plot.width = 8,repr.plot.height =  8)
rownames(gsva.kegg) = gsub('hallmark_','',rownames(gsva.kegg))


l <- data.frame()
cells = unique(scRNA$Celltype)
for (i in cells){
  Idents(scRNA) = 'Celltype'
  cell = subset(scRNA,idents = i)
  DefaultAssay(cell) = 'RNA'
  Idents(cell) = 'Group'
  diffgene = FindMarkers(cell,ident.1 = 'Tumor',logfc.threshold = 0.25,
                         test.use = "wilcox",
                         min.pct = 0.15)
  diffgene$cluster = i
  diffgene$gene = rownames(diffgene)
  l <- rbind(l, diffgene)
}

head(l,2)

sce.markers = l %>%  filter(avg_log2FC > 0.15)

ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster)

xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
#如果原始的ID号为entrez gene id那么这里keyType设置为ENTREZID
xx<-setReadable(xx, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# down
sce.markers = l %>%  filter(avg_log2FC < -0.15)

ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster)

xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)

#如果原始的ID号为entrez gene id那么这里keyType设置为ENTREZID
xx<-setReadable(xx, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


# All
for (i in unique(scRNA$Celltype) ){

  Idents(scedata) = 'Celltype'
  cell = subset(scedata,Celltype == i)
  
  DefaultAssay(cell) = 'RNA'
  Idents(cell) = 'Group'
  diffgene = FindMarkers(cell,ident.1 = 'Tumor' ) %>%  filter(p_val<0.05 )
  diffgene$gene  = rownames(diffgene)

  ego_CC <- enrichGO(gene          = diffgene$gene ,
                     OrgDb         = 'org.Hs.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.5)
  ego_cc <- data.frame(ego_CC)
  ego_MF <- enrichGO(gene          = diffgene$gene ,
                     OrgDb         = 'org.Hs.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.5)
  ego_mf <- data.frame(ego_MF)
  ego_BP <- enrichGO(gene          = diffgene$gene ,
                     OrgDb         = 'org.Hs.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.5) 
  ego_bp <- data.frame(ego_BP)
  
  
  genelist <- bitr(diffgene$gene, fromType="SYMBOL",
                   toType="ENTREZID", OrgDb='org.Hs.eg.db')
  genelist <- pull(genelist,ENTREZID)               
  ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
  
}



