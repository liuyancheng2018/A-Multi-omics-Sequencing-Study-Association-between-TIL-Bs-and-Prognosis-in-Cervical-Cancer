library(DoubletFinder)
library(dplyr)
library(Seurat)
library(ggplot2)


sample = snakemake@params[[1]]

run_Finddoublet <- function(sample_name, n_capture_cell){
  n_capture_cell = as.character(n_capture_cell)
  multiplet_rate = c('500' = 0.004, '1000' = 0.008, '2000' = 0.016, '3000' = 0.023, '4000' = 0.031,
                     '5000' = 0.039, '6000' = 0.046, '7000' = 0.054, '8000' = 0.061, '9000' = 0.069,
                     '10000' = 0.076)
  stopifnot(n_capture_cell %in% names(multiplet_rate))
  single_data = Read10X(paste0("data/", sample_name, "/filtered_feature_bc_matrix/"))
  seu_obj = CreateSeuratObject(counts = single_data, project = sample_name)
  seu_obj = SCTransform(seu_obj)
  seu_obj <- RunPCA(seu_obj, verbose = F)
  seu_obj <- RunUMAP(seu_obj, dims = 1:15)
  seu_obj <- FindNeighbors(seu_obj, dims = 1:15)
  seu_obj <- FindClusters(seu_obj, resolution = 0.5)
  n_cpu = 1
  sweep.res.list = paramSweep_v3(seu_obj, PCs = 1:15, sct = TRUE, num.cores = n_cpu)
  sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pn = 0.25
  pk = bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  nexp = round(multiplet_rate[n_capture_cell]*nrow(seu_obj@meta.data)) %>% unname()
  seu_obj = doubletFinder_v3(seu_obj, PCs = 1:15, pN = pn, pK = pk, nExp = nexp, reuse.pANN = FALSE, sct = TRUE)
  DimPlot(seu_obj, reduction = 'umap', group.by = paste0("DF.classifications_", pn, "_", pk, "_", nexp))
  ggsave(filename = paste0("doubltremove/DimPlot_",sample_name,"_doublet_cell.pdf"))
  data = seu_obj[[paste0("DF.classifications_", pn, "_", pk, "_", nexp)]]
  doublet_barcode = as.data.frame(rownames(filter(data, data[,1] == "Doublet")))
  order = which(sample == sample_name)
  doublet_barcode$tmp = paste(doublet_barcode[, 1], "_", order, sep = "")
  write.table(doublet_barcode[, 2], 
              paste0("doubltremove/", sample_name, "_doublet_barcode.txt"), 
              sep = "\t", row.names = F, col.names = F)
  return
}

for (i in sample){
  run_Finddoublet(sample_name = i, n_capture_cell = 10000)
}
