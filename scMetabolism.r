################################################
################################################
### Aims:
### 1. scMetabolism analysis of CDH11 epithelial cells
### 2. scMetabolism analysis of ITGA3 epithelial cells
### 3. scMetabolism analysis of PTGDS epithelial cells

######## 1. scMetabolism analysis of CDH11 epithelial cells ########
library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(ggsci)
library(tidyverse)
library(Seurat)
epi <- readRDS("./scRNA_epi_1009_noprcc.rds")
table(epi$Group)
table(epi$Celltype3)
table(epi$new_sample)
epi@meta.data$Celltype3 <- recode(epi@meta.data$Celltype3,
                                  "Type_B_Intercalated_cell" = "ICB")
scobj <- epi
Idents(scobj) <- "Celltype3"
scobj <- subset(scobj, idents = c("CDH11_epi"))
scobj$Celltype3 <- as.character(scobj$Celltype3)
table(scobj$Group, scobj$Celltype3)
Idents(scobj) <- "orig.ident"
countexp.Seurat <- sc.metabolism.Seurat(obj = scobj,  
                                        method = "AUCell", 
                                        imputation = F, 
                                        ncores = 10, 
                                        metabolism.type = "KEGG")
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change))
length(rownames(countexp.Seurat@assays$METABOLISM$score))
rownames(score) = gsub(' / ','_',rownames(score) )
colnames(countexp.Seurat@meta.data) = gsub(' / ','_',colnames(countexp.Seurat@meta.data) )
library(ggpubr)
dim(countexp.Seurat@meta.data)
MY =list(c('Tumor','Adjacent'))
Idents(countexp.Seurat) <- "Group"


######## 2. scMetabolism analysis of ITGA3 epithelial cells ########
epi <- readRDS("./scRNA_epi_1009_noprcc.rds")
table(epi$Group)
table(epi$Celltype3)
table(epi$new_sample)
epi@meta.data$Celltype3 <- recode(epi@meta.data$Celltype3,
                                  "Type_B_Intercalated_cell" = "ICB")
scobj <- epi
Idents(scobj) <- "Celltype3"
scobj <- subset(scobj, idents = c("ITGA3_epi"))
scobj$Celltype3 <- as.character(scobj$Celltype3)
table(scobj$Group, scobj$Celltype3)
Idents(scobj) <- "orig.ident"
countexp.Seurat <- sc.metabolism.Seurat(obj = scobj,
                                        method = "AUCell", 
                                        imputation = F, 
                                        ncores = 10, 
                                        metabolism.type = "KEGG")
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change))
length(rownames(countexp.Seurat@assays$METABOLISM$score))
rownames(score) = gsub(' / ','_',rownames(score) )
colnames(countexp.Seurat@meta.data) = gsub(' / ','_',colnames(countexp.Seurat@meta.data) )
library(ggpubr)
dim(countexp.Seurat@meta.data)
MY =list(c('Tumor','Adjacent'))
Idents(countexp.Seurat) <- "Group"


######## 3. scMetabolism analysis of PTGDS epithelial cells ########
epi <- readRDS("./scRNA_epi_1009_noprcc.rds")
table(epi$Group)
table(epi$Celltype3)
table(epi$new_sample)
epi@meta.data$Celltype3 <- recode(epi@meta.data$Celltype3,
                                  "Type_B_Intercalated_cell" = "ICB")
scobj <- epi
Idents(scobj) <- "Celltype3"
scobj <- subset(scobj, idents = c("PTGDS_epi"))
scobj$Celltype3 <- as.character(scobj$Celltype3)
table(scobj$Group, scobj$Celltype3)
Idents(scobj) <- "orig.ident"
countexp.Seurat <- sc.metabolism.Seurat(obj = scobj,
                                        method = "AUCell", 
                                        imputation = F, 
                                        ncores = 10, 
                                        metabolism.type = "KEGG")
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change))
length(rownames(countexp.Seurat@assays$METABOLISM$score))
rownames(score) = gsub(' / ','_',rownames(score) )
colnames(countexp.Seurat@meta.data) = gsub(' / ','_',colnames(countexp.Seurat@meta.data) )
library(ggpubr)
dim(countexp.Seurat@meta.data)
MY =list(c('Tumor','Adjacent'))
Idents(countexp.Seurat) <- "Group"

