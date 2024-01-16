library(Seurat)
library(SeuratObject)
library(Scissor)
library(future)
library(preprocessCore)

options(future.globals.maxSize = 10 * 1024^3)
load(snakemake@input[[1]])

DefaultAssay(sce.SCTnormalize.clean) <- "RNA"
data = sce.SCTnormalize.clean[["RNA"]]@counts
#data <- as.matrix(data)
sc_dataset <- Seurat_preprocessing(data, verbose = F)

bulk_dataset = read.table(snakemake@input[[2]], sep = "\t", header = T, row.names = 1)
phenotype = read.table(snakemake@input[[3]], sep = "\t", header = T, row.names = 1)

survival = phenotype[, c("days_to_last_follow_up", "vital_status")]
rownames(survival) <- gsub("-", ".", rownames(survival))
colnames(survival) <- c("time", "status")
survival <- survival[colnames(bulk_dataset), ]
rownames(survival) <- NULL
plan("multicore", workers = 24)
infos1 <- Scissor(as.matrix(bulk_dataset), 
                  sc_dataset, 
                  survival, 
                  alpha = 0.05, 
                  family = "cox", 
                  Save_file = snakemake@output[[1]])

save(sce.SCTnormalize.clean, bulk_dataset, sc_dataset, phenotype, survival, infos1, file = snakemake@output[[2]])
