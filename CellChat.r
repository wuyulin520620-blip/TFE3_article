################################################
################################################
### Aims:
### 1. CellChat (All Cell Types)
### 2. CellChat (subtypes of T cells, epithelial cells, and macrophages)
### 3. CellChat (subtypes of epithelial cells and fibroblasts)

######## 1. CellChat (All Cell Types) ########
library(CellChat)
library(patchwork)
sce <- readRDS("./scRNA_merge_0930_noprcc.rds")
Idents(sce) <- "Celltype"
table(sce$Celltype)

table(sce$Group)
wt <- subset(sce, Group == "Adjacent")
ko <- subset(sce, Group == "Tumor")

wt <- createCellChat(wt@assays$RNA@data, meta = wt@meta.data, group.by = "Celltype")
ko <- createCellChat(ko@assays$RNA@data, meta = ko@meta.data, group.by = "Celltype")

length(unique(sce$Celltype))

cellchat <- wt
cellchat@DB <-  CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
wt <- cellchat
unique(ko@idents)

cellchat <- ko
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ko <- cellchat

saveRDS(ko, file = paste0(output, "ko_cellchat_all.rds"))
saveRDS(wt, file = paste0(output, "wt_cellchat_all.rds"))


######## 2. CellChat (subtypes of T cells, epithelial cells, and macrophages) ########
library(CellChat)
library(patchwork)
library(Seurat)
sce <- readRDS("./scRNA_merge_0930_noprcc.rds")
epi <- readRDS("./scRNA_epi_1009_noprcc.rds")
table(epi$Group)
table(epi$Celltype3)
table(epi$new_sample)
epi@meta.data$Celltype3 <- recode(epi@meta.data$Celltype3,
                                  "Type_B_Intercalated_cell" = "ICB")
t <- readRDS("./sub.T.final.rds")
table(t$group)
table(t$Celltype)
table(t$orig.ident)
macro <- readRDS("./scRNA_mye_1013.rds")
table(macro$Group)
table(macro$Celltype2)
table(macro$new_sample)

# T_annotaion
t$Celltype <- as.character(t$Celltype)
match_idx <- match(rownames(sce@meta.data), rownames(t@meta.data))
has_match <- !is.na(match_idx)
sce@meta.data$Celltype1[has_match] <- t@meta.data[match_idx[has_match], "Celltype"]

# Epi_annotaion
epi$Celltype3 <- as.character(epi$Celltype3)
match_idx <- match(rownames(sce@meta.data), rownames(epi@meta.data))
has_match <- !is.na(match_idx)
sce@meta.data$Celltype1[has_match] <- epi@meta.data[match_idx[has_match], "Celltype3"]

# Macrophage_annotaion
macro$Celltype2 <- as.character(macro$Celltype2)
match_idx <- match(rownames(sce@meta.data), rownames(macro@meta.data))
has_match <- !is.na(match_idx)
sce@meta.data$Celltype1[has_match] <- macro@meta.data[match_idx[has_match], "Celltype2"]

table(sce@meta.data$Celltype1)
celltypes <- unique(sce@meta.data$Celltype1)
cat('celltypes <- c("', paste(celltypes, collapse = '", "'), '")\n', sep = "")
Idents(sce) <- "Celltype1"

table(sce@meta.data$Celltype)

table(sce@meta.data$new_sample)
celltypes <- unique(sce@meta.data$new_sample)
cat('celltypes <- c("', paste(celltypes, collapse = '", "'), '")\n', sep = "")
table(sce@meta.data$Group)
Idents(sce) <- "Group"
Adjacent <- subset(sce, idents = c("Adjacent"))

Idents(sce) <- "new_sample"
sce <- subset(sce, idents = c("ASPL_T2", "ASPL_T1", "MED15_T1", 
                              "MED15_T2", "NONO_T1", "2T", "3T"))
sce@meta.data$new_sample <- recode(sce@meta.data$new_sample,
                                   "2T" = "NONO_2T",
                                   "3T" = "ASPL_3T")
sce@meta.data$Group1 <- recode(sce@meta.data$new_sample,
                               "ASPL_3T" = "ASPL",
                               "ASPL_T1" = "ASPL",
                               "ASPL_T2" = "ASPL",
                               "MED15_T1" = "MED15",
                               "MED15_T2" = "MED15",
                               "NONO_2T" = "NONO",
                               "NONO_T1" = "NONO"
)
sce@meta.data$Group1 <- as.character(sce@meta.data$Group1)
table(sce@meta.data$Group1)
saveRDS(sce, file = paste0(output, "scRNA_ASPL_MED15_NONO.rds"))

table(sce$Group1)
ASPL <- subset(sce, Group1 == "ASPL")
MED15 <- subset(sce, Group1 == "MED15")
MED15$Celltype1 <- as.character(MED15$Celltype1)
NONO <- subset(sce, Group1 == "NONO")

Adjacent <- createCellChat(Adjacent@assays$RNA@data, meta = Adjacent@meta.data, group.by = "Celltype")
ASPL <- createCellChat(ASPL@assays$RNA@data, meta = ASPL@meta.data, group.by = "Celltype")
MED15 <- createCellChat(MED15@assays$RNA@data, meta = MED15@meta.data, group.by = "Celltype")
NONO <- createCellChat(NONO@assays$RNA@data, meta = NONO@meta.data, group.by = "Celltype")

length(unique(sce$Celltype))

cellchat <- Adjacent
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
Adjacent <- cellchat
unique(Adjacent@idents)

cellchat <- ASPL
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ASPL <- cellchat
unique(ASPL@idents)

cellchat <- MED15
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
MED15 <- cellchat

cellchat <- NONO
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
NONO <- cellchat

saveRDS(Adjacent, file = paste0(output, "Adjacent_cellchat_all.rds"))
saveRDS(ASPL, file = paste0(output, "ASPL_cellchat_all.rds"))
saveRDS(MED15, file = paste0(output, "MED15_cellchat_all.rds"))
saveRDS(NONO, file = paste0(output, "NONO_cellchat_all.rds"))

Adjacent <- readRDS(file = paste0(output, "Adjacent_cellchat_all.rds"))
ASPL <- readRDS(file = paste0(output, "ASPL_cellchat_all.rds"))
MED15 <- readRDS(file = paste0(output, "MED15_cellchat_all.rds"))
NONO <- readRDS(file = paste0(output, "NONO_cellchat_all.rds"))

object.list <- list(Adjacent = Adjacent, ASPL = ASPL, MED15 = MED15 , NONO = NONO)

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)


######## 3. CellChat (subtypes of epithelial cells and fibroblasts) ########
library(CellChat)
library(patchwork)
library(Seurat)
library(qs)
sce <- readRDS("./scRNA_merge_0930_noprcc.rds")
sce$Celltype1 <- sce$Celltype
sce$Celltype1 <- as.character(sce$Celltype1)
fib <- readRDS("./scRNA_Fib_1010_noprcc.rds")
table(fib$Celltype)
epi <- readRDS("./scRNA_epi_1009_noprcc.rds")
table(epi$Group)
table(epi$Celltype3)
table(epi$new_sample)
epi@meta.data$Celltype3 <- recode(epi@meta.data$Celltype3,
                                  "Type_B_Intercalated_cell" = "ICB")
# Fib_annotaion
fib$Celltype <- as.character(fib$Celltype)
match_idx <- match(rownames(sce@meta.data), rownames(fib@meta.data))
has_match <- !is.na(match_idx)
sce@meta.data$Celltype1[has_match] <- fib@meta.data[match_idx[has_match], "Celltype"]
table(sce@meta.data$Celltype1)
# Epi_annotaion
epi$Celltype3 <- as.character(epi$Celltype3)
match_idx <- match(rownames(sce@meta.data), rownames(epi@meta.data))
has_match <- !is.na(match_idx)
sce@meta.data$Celltype1[has_match] <- epi@meta.data[match_idx[has_match], "Celltype3"]

table(sce@meta.data$Celltype1)
celltypes <- unique(sce@meta.data$Celltype1)
cat('celltypes <- c("', paste(celltypes, collapse = '", "'), '")\n', sep = "")
Idents(sce) <- "Celltype1"
sce <- subset(sce, idents = c("SPP1_epi", "ITGA3_epi", "CDH11_epi", "PTGDS_epi", "ICB", "TAL", "PT", "PC", 
                              "Collagen_Fib", "NACM1_Fib"))
sce@meta.data$Celltype1 <- factor(sce@meta.data$Celltype1, levels = c("SPP1_epi", "ITGA3_epi", "CDH11_epi", "PTGDS_epi", "ICB", "TAL", "PT", "PC", 
                                                                      "Collagen_Fib", "NACM1_Fib"))
table(sce@meta.data$Celltype1)
table(sce@meta.data$new_sample)
Idents(sce) <- "Group"

table(sce$Group)
Idents(sce) <- "Group"
Adjacent <- subset(sce, idents = c("Adjacent"))

Idents(sce) <- "new_sample"
sce <- subset(sce, idents = c("ASPL_T2", "ASPL_T1", "MED15_T1", 
                              "MED15_T2", "NONO_T1", "2T", "3T"))
sce@meta.data$new_sample <- recode(sce@meta.data$new_sample,
                                   "2T" = "NONO_2T",
                                   "3T" = "ASPL_3T")
sce@meta.data$Group1 <- recode(sce@meta.data$new_sample,
                               "ASPL_3T" = "ASPL",
                               "ASPL_T1" = "ASPL",
                               "ASPL_T2" = "ASPL",
                               "MED15_T1" = "MED15",
                               "MED15_T2" = "MED15",
                               "NONO_2T" = "NONO",
                               "NONO_T1" = "NONO"
)
sce@meta.data$Group1 <- as.character(sce@meta.data$Group1)
table(sce@meta.data$Group1)

table(sce$Group1)
ASPL <- subset(sce, Group1 == "ASPL")
MED15 <- subset(sce, Group1 == "MED15")
MED15$Celltype1 <- as.character(MED15$Celltype1)
NONO <- subset(sce, Group1 == "NONO")

Adjacent <- createCellChat(Adjacent@assays$RNA@data, meta = Adjacent@meta.data, group.by = "Celltype1")
ASPL <- createCellChat(ASPL@assays$RNA@data, meta = ASPL@meta.data, group.by = "Celltype1")
MED15 <- createCellChat(MED15@assays$RNA@data, meta = MED15@meta.data, group.by = "Celltype1")
NONO <- createCellChat(NONO@assays$RNA@data, meta = NONO@meta.data, group.by = "Celltype1")

length(unique(sce$Celltype1))

cellchat <- Adjacent
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
Adjacent <- cellchat
unique(Adjacent@idents)

cellchat <- ASPL
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ASPL <- cellchat
unique(ASPL@idents)

cellchat <- MED15
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
MED15 <- cellchat

cellchat <- NONO
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
NONO <- cellchat

saveRDS(Adjacent, file = paste0(output_dir, "Adjacent_cellchat_all.rds"))
saveRDS(ASPL, file = paste0(output_dir, "ASPL_cellchat_all.rds"))
saveRDS(MED15, file = paste0(output_dir, "MED15_cellchat_all.rds"))
saveRDS(NONO, file = paste0(output_dir, "NONO_cellchat_all.rds"))

#Adjacent <- readRDS(file = paste0(output, "Adjacent_cellchat_all.rds"))
ASPL <- readRDS(file = paste0(output_dir, "ASPL_cellchat_all.rds"))
MED15 <- readRDS(file = paste0(output_dir, "MED15_cellchat_all.rds"))
NONO <- readRDS(file = paste0(output_dir, "NONO_cellchat_all.rds"))

object.list <- list(ASPL = ASPL, MED15 = MED15 , NONO = NONO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)