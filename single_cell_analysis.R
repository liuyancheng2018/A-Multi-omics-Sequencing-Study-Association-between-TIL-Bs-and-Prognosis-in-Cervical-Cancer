set.seed(1234)

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

load_Rdata = snakemake@params[[1]]
cluster2celltype = snakemake@params[[2]]

color_file = read.table(snakemake@input[[1]], sep = "\t", comment.char = "", header = T)
use_color = color_file$color
names(use_color) = color_file$cell_type

load(load_Rdata)
print(cluster2celltype)
sce.SCTnormalize.clean = sce.SCTnormalize.harmony
sce.SCTnormalize.clean[['cell_type']] = unname(cluster2celltype[sce.SCTnormalize.clean@meta.data$seurat_clusters]) 
colnames(sce.SCTnormalize.clean@meta.data)


pdf(file = snakemake@output[[1]])
DimPlot(sce.SCTnormalize.clean, reduction = 'umap', group.by = 'cell_type',
        label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

pdf(file = snakemake@output[[2]])
DimPlot(sce.SCTnormalize.clean, reduction = 'tsne', group.by = 'cell_type',
        label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

pdf(file = snakemake@output[[3]])
DimPlot(sce.SCTnormalize.clean, reduction = 'umap', group.by = 'orig.ident',
        label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

pdf(file = snakemake@output[[4]])
DimPlot(sce.SCTnormalize.clean, reduction = 'tsne', group.by = 'orig.ident',
        label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

Idents(sce.SCTnormalize.clean) = "cell_type"
tsne_data = DimPlot(sce.SCTnormalize.clean, reduction = "tsne")
tsne_data = tsne_data$data

p1 = ggplot()+geom_point(data = tsne_data, mapping = aes(tSNE_1, tSNE_2, color = ident), size = 0.5) +
  scale_color_manual(values = use_color) +
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

ggsave(snakemake@output[[11]], width = 8, height = 8) ##

p2 = p1 + geom_segment(aes(x = min(tsne_data$tSNE_1), y = min(tsne_data$tSNE_2),
                           xend = min(tsne_data$tSNE_1) + 15, yend = min(tsne_data$tSNE_2)),
                       colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(tsne_data$tSNE_1), y = min(tsne_data$tSNE_2),
                   xend = min(tsne_data$tSNE_1), yend = min(tsne_data$tSNE_2) + 15),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(tsne_data$tSNE_1) + 7, y = min(tsne_data$tSNE_2) - 2, label = "tSNE_1",
           color = "black", size = 5, fontface = "bold" ) + 
  annotate("text", x = min(tsne_data$tSNE_1) - 2, y = min(tsne_data$tSNE_2) + 7, label = "tSNE_2",
           color = "black", size = 5, fontface = "bold" ,angle=90)

ggsave(snakemake@output[[5]], width = 8, height = 8)

p3 = p1 + geom_segment(aes(x = min(tsne_data$tSNE_1), y = min(tsne_data$tSNE_2),
                           xend = min(tsne_data$tSNE_1) + 15, yend = min(tsne_data$tSNE_2)),
                       colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(tsne_data$tSNE_1), y = min(tsne_data$tSNE_2),
                   xend = min(tsne_data$tSNE_1), yend = min(tsne_data$tSNE_2) + 15),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))

ggsave(snakemake@output[[6]], width = 8, height = 8)

Idents(sce.SCTnormalize.clean) = "cell_type"
umap_data = DimPlot(sce.SCTnormalize.clean, reduction = "umap")
umap_data = umap_data$data

p4 = ggplot()+geom_point(data = umap_data, mapping = aes(UMAP_1, UMAP_2, color = ident), size = 0.5) +
  scale_color_manual(values = use_color) +
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position="none")

ggsave(snakemake@output[[10]], width = 8, height = 8) ##

p5 = p4 + geom_segment(aes(x = min(umap_data$UMAP_1), y = min(umap_data$UMAP_2),
                           xend = min(umap_data$UMAP_1) + 4, yend = min(umap_data$UMAP_2)),
                       colour = "black", size = 1, arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap_data$UMAP_1), y = min(umap_data$UMAP_2),
                   xend = min(umap_data$UMAP_1), yend = min(umap_data$UMAP_2) + 5),
               colour = "black", size = 1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap_data$UMAP_1) + 2, y = min(umap_data$UMAP_2) - 1, label = "UMAP_1",
           color = "black", size = 8, fontface = "bold" ) + 
  annotate("text", x = min(umap_data$UMAP_1) - 1, y = min(umap_data$UMAP_2) + 2, label = "UMAP_2",
           color = "black", size = 8, fontface = "bold" ,angle=90)

ggsave(snakemake@output[[7]], width = 8, height = 8)

p6 = p4 + geom_segment(aes(x = min(umap_data$UMAP_1), y = min(umap_data$UMAP_2),
                           xend = min(umap_data$UMAP_1) + 4, yend = min(umap_data$UMAP_2)),
                       colour = "black", size = 1, arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap_data$UMAP_1), y = min(umap_data$UMAP_2),
                   xend = min(umap_data$UMAP_1), yend = min(umap_data$UMAP_2) + 5),
               colour = "black", size = 1, arrow = arrow(length = unit(0.3,"cm")))

ggsave(snakemake@output[[8]], width = 8, height = 8)

save(sce.SCTnormalize.clean, file = snakemake@output[[9]])

file_path = snakemake@params[[3]]

for (i in 1:nrow(color_file)){
  j <- color_file$cell_type[i]
  color <- color_file$color[i]
  pdf(file = paste0(file_path, "/label_", j, ".pdf"), width = 7, height = 7)
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, xlim = c(-1, 1), ylim = c(-1, 1), type = "n", asp = 1,
       xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)
  center <- c(0, 0)
  radius <- 1
  symbols(center[1], center[2], circles = radius,
          inches = FALSE, add = TRUE, bg = color, lty = 0)
  dev.off()
}

