################################################
################################################
### Aims:
### 1. monocle2

######## 1. monocle2 ########
library(progeny)
library(Seurat)
library(monocle)
library(tidyverse)
epi <- readRDS(paste0(input_dir, "scRNA_epi_1009_noprcc.rds"))
table(epi$Group)
table(epi$Celltype3)
table(epi$new_sample)
epi@meta.data$Celltype3 <- recode(epi@meta.data$Celltype3,
                                  "Type_B_Intercalated_cell" = "ICB")
epi$new_sample <- as.character(epi$new_sample)
scobj <- subset(epi, Group == "Tumor")

table(scobj$new_sample)
pathways <- progeny(as.matrix(scobj@assays$RNA@data), scale = TRUE, organism = "Human", top = 100, perm = 1)
pathways <- as.data.frame(pathways)
head(pathways)
scobj <- scobj[,rownames(pathways)]
scobj <- AddMetaData(scobj, metadata = pathways)

scobj@meta.data$new_sample <- recode(scobj@meta.data$new_sample,
                                     "2T" = "NONO_2T",
                                     "3T" = "ASPL_3T")
scobj@meta.data$Group1 <- recode(scobj@meta.data$new_sample,
                                 "ASPL_3T" = "ASPL",
                                 "ASPL_T1" = "ASPL",
                                 "ASPL_T2" = "ASPL",
                                 "MED15_T1" = "MED15",
                                 "MED15_T2" = "MED15",
                                 "NONO_2T" = "NONO",
                                 "NONO_T1" = "NONO")

expr_matrix <- as(as.matrix(scobj@assays$RNA@counts), 'sparseMatrix')
p_data <- scobj@meta.data 
f_data <- data.frame(gene_short_name = row.names(scobj),row.names = row.names(scobj))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size()) 
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds)
saveRDS(cds, paste0(output_dir, "cds.rds"))
cds <- readRDS(paste0(output_dir, "cds.rds"))
cds <- orderCells(cds, root_state = 2)
plot_cell_trajectory(cds)
