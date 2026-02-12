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
library(DoubletFinder)

color1 = ggsci::pal_nejm(palette = c("default"), alpha = 0.8)(8)
color2 = ggsci::pal_aaas(palette = c("default"), alpha = 0.8)(8)
color3 = ggsci::pal_jama(palette = c("default"), alpha = 0.8)(8)

color = c(color1,color2,color3)

matrix_to_seurat <- function(file_id){
  file_name = paste("./01.rawdata/",file_id,"_matrix.tsv.gz",sep ="")
  counts=fread(file = file_name,data.table = T,sep = '\t')
  counts=data.frame(counts)
  rownames(counts)=counts[,1]
  counts=counts[,-1]
  colnames(counts)=paste0(file_id,"_",colnames(counts)) # 在细胞前添加样本名 
  genename <- rownames(counts)
  genename<-gsub("_","-",genename) # 将基因名中"_"替换为"-",防止seurat识别的更改
  rownames(counts)=genename
  sce <- CreateSeuratObject(counts=counts,project = file_id,min.cells = 3, min.features = 250) # 构建seurat对象，保留至少在3个细胞中表达的基因 & 至少有250基因表达的细胞
  return(sce)
}

#### 读数据 ####
if (platform=="10x") {
  folders <- list.files("./01.rawdata/")
  datalist=list()
  for (i in 1:length(folders)){
    datalist[[i]]<- CreateSeuratObject(counts = Read10X(data.dir = paste0("./01.rawdata/",folders[i])), 
                                       project = folders[i], 
                                       min.cells = 3, 
                                       min.features = 250)}
  sce <- merge(datalist[[1]],y=datalist[2:length(datalist)],add.cell.ids = folders)
  datalist <- SplitObject(sce,split.by = "orig.ident")
  names(datalist)=folders
  rm(sce)
  
} else {
  file_id <- list.files("./01.rawdata/",pattern='matrix.tsv.gz$', full=FALSE)
  file_id<- gsub("_matrix.tsv.gz","",file_id)
  datalist=list()
  for (i in 1:length(file_id)){
    datalist[[i]]=matrix_to_seurat(file_id[i])}
  names(datalist)=file_id
}

for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比
  datalist[[i]] <- sce
  rm(sce)
}

sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
sce <- subset(sce, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

#### 合并降维聚类 ##########
sce <-NormalizeData(sce)
sce <- FindVariableFeatures(object = sce,selection.method = 'vst', nfeatures = 3000)
ifnb.list <- SplitObject(sce, split.by = "orig.ident")
seu.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:50)
seu.integrated <- IntegrateData(anchorset = seu.anchors, dims = 1:50)
DefaultAssay(seu.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
seu.integrated <- ScaleData(seu.integrated, features = rownames(seu.integrated))
seu.integrated <- RunPCA(seu.integrated, npcs = 50, verbose = FALSE)
seu.integrated <- RunUMAP(seu.integrated, dims = 1:30)
seu.integrated <- RunTSNE(seu.integrated, dims = 1:30)
seu.integrated <- FindNeighbors(seu.integrated, dims = 1:30)
seu.integrated <- FindClusters(seu.integrated, resolution = 0.2)

#### 去双细胞 ##########
sweep.res.list_kidney <- paramSweep_v3(seu.integrated, PCs = 1:20, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
pK_bcmvn <- bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(seu.integrated)*8*1e-7

homotypic.prop <- modelHomotypic(seu.integrated$orig.ident)
nExp_poi <- round(DoubletRate*ncol(seu.integrated)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

doubletFinder = function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
                          sct = FALSE, annotations = NULL) 
{
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", 
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    orig.commands <- seu@commands
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                                 model.use = orig.commands$ScaleData.RNA$model.use, 
                                 do.scale = orig.commands$ScaleData.RNA$do.scale, 
                                 do.center = orig.commands$ScaleData.RNA$do.center, 
                                 scale.max = orig.commands$ScaleData.RNA$scale.max, 
                                 block.size = orig.commands$ScaleData.RNA$block.size, 
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                              npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                              verbose = FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                 ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                             ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbors_that_are_doublets = neighbors[neighbors > 
                                                    n_real.cells]
          if (length(neighbors_that_are_doublets) > 0) {
            neighbor_types[i, ] <- table(doublet_types1[neighbors_that_are_doublets - 
                                                          n_real.cells]) + table(doublet_types2[neighbors_that_are_doublets - 
                                                                                                  n_real.cells])
            neighbor_types[i, ] <- neighbor_types[i, 
            ]/sum(neighbor_types[i, ])
          }
          else {
            neighbor_types[i, ] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 
                                                                    1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) = levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors", 
                              pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, 
                                                                              ct]
      }
    }
    return(seu)
  }
}

scRNA <- doubletFinder(seu.integrated, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, 
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
scRNA = subset(scRNA, subset = DF.classifications_0.25_0.005_2857  != 'Doublet')

scRNA = RenameIdents(scRNA,
                     '0' = 'Epithelial cell',
                     '1' = 'Epithelial cell',
                     '2' = 'T cell',
                     '3' = 'Macrophage',
                     '4' = 'Fibroblast',
                     '5' = 'Endothelial cell',
                     '6' = 'Epithelial cell',
                     '7' = 'Epithelial cell',
                     '8' = 'B cell',
                     '9' = 'Epithelial cell',
                     '10' = 'Podocyte',
                     '11' = 'Prof',
                     '12' = 'Plasma',
                     '13' = 'Neutrophils',
                     '14' = 'Mast',
                     '15' = 'Endothelial cell')
scRNA <- AddMetaData(scRNA,scRNA@active.ident,col.name = "Celltype")
saveRDS(scRNA,'scRNA.rds')