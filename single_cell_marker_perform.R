library(Seurat)
library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(tidyverse)

Load_rdata = snakemake@params[[1]]
Order = snakemake@params[[3]]
Figure_path = snakemake@params[[4]]

load(Load_rdata)

markergenes = read.table(snakemake@input[[1]], sep = "\t", header = T)

Idents(sce.SCTnormalize.clean) = "cell_type"

tmp_data <- melt(markergenes, id.vars = c("cluster"), variable.name = "marker", value.name = "gene")
bubble_data1 <- ordered(tmp_data$cluster, levels=Order)
bubble_data2 = tmp_data[order(bubble_data1),][, c(1, 3)]

unique(sce.SCTnormalize.clean$cell_type)
length(unique(sce.SCTnormalize.clean$cell_type))
sce.SCTnormalize.clean$cell_type = factor(sce.SCTnormalize.clean$cell_type ,levels = Order)

pdf(file = snakemake@output[[1]], width = 2.5, height = 4)
print(DotPlot(sce.SCTnormalize.clean, features = bubble_data2$gene) + RotatedAxis() + coord_flip() +
        scale_x_discrete("") + scale_y_discrete(limits = Order) +
        scale_color_gradient2(low = 'navy', high = 'firebrick3', mid = "white"))
dev.off()

pdf(file = snakemake@output[[2]], width = 2.5, height = 4)
DotPlot(sce.SCTnormalize.clean, features = bubble_data2$gene) + RotatedAxis() + coord_flip() + 
  scale_x_discrete("") + scale_y_discrete(limits = Order) + xlab("") + ylab("") +
  scale_color_gradient2(low = 'navy', high = 'firebrick3', mid = "white") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")
dev.off()

pdf(file = snakemake@output[[3]], width = 12, height = 5)
DotPlot(sce.SCTnormalize.clean, features = bubble_data2$gene) + RotatedAxis() + 
        scale_x_discrete("") + scale_y_discrete(limits = Order) + 
        scale_color_gradient2(low = 'navy', high = 'firebrick3', mid = "white")
dev.off()

pdf(file = snakemake@output[[4]], width = 12, height = 5)
DotPlot(sce.SCTnormalize.clean, features = bubble_data2$gene) + RotatedAxis() + 
  scale_x_discrete("") + scale_y_discrete(limits = Order) + xlab("") + ylab("") +
  scale_color_gradient2(low = 'navy', high = 'firebrick3', mid = "white") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")
dev.off()

pdf(file = snakemake@output[[5]])
VlnPlot(sce.SCTnormalize.clean, features = bubble_data2$gene, stack = TRUE) &
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(angle = 45, size = 10, hjust = 0, vjust = 0),
        legend.position = "none")
dev.off()

pdf(file = snakemake@output[[6]])
VlnPlot(sce.SCTnormalize.clean, features = bubble_data2$gene, stack = TRUE) &
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x.top = element_blank(),
        legend.position = "none")
dev.off()

pdf(file = snakemake@output[[7]])
VlnPlot(sce.SCTnormalize.clean, features = bubble_data2$gene, stack = TRUE, flip = TRUE) + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

pdf(file = snakemake@output[[8]])
VlnPlot(sce.SCTnormalize.clean, features = bubble_data2$gene, stack = TRUE, flip = TRUE) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text.y.left = element_blank(),
        strip.text.y.right = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()) + xlab("") + ylab("")
dev.off()

features = bubble_data2$gene
Norm_data <- FetchData(sce.SCTnormalize.clean, features, slot = "data")

Norm_data$Cell <- colnames(sce.SCTnormalize.clean)
Norm_data$Idents <- sce.SCTnormalize.clean@meta.data$cell_type

Norm_data <- reshape2::melt(Norm_data, id.vars = c("Cell","Idents"), measure.vars = features,
                            variable.name = "Feat", value.name = "Expr")


pdf(file = snakemake@output[[9]])
ggplot(Norm_data, aes(factor(Idents), Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5)) + xlab("") + ylab("")
dev.off()


pdf(file = snakemake@output[[10]])
ggplot(Norm_data, aes(factor(Idents), Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text.y.left = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()) + xlab("") + ylab("")
dev.off()


markergenes1 = read.table(snakemake@input[[1]], sep = "\t", header = T, row.names = 1)


for (i in Order){
  markergenes2 = markergenes1[i, ] 
  markergenes2 = as.character(markergenes2)
  Norm_data <- FetchData(sce.SCTnormalize.clean, markergenes2, slot = "data")
  Norm_data$means = rowMeans(Norm_data)
  umap_plot = sce.SCTnormalize.clean@reductions$umap@cell.embeddings %>% as.data.frame()
  new_data = cbind(umap_plot, Norm_data$means)
  title_tmp = paste0("(", paste0(markergenes2, collapse = ", "),")")
  pdf(paste0(Figure_path, "/", i, "_marker_average_expression_UMAP.pdf"))
  print(ggplot(new_data, aes(x = UMAP_1, y = UMAP_2, colour = Norm_data$means)) + 
          geom_point(size = 0.2, alpha = 1) + 
          scale_colour_gradientn(colours = c("#CCCCCC", "red")) + 
          labs(title = paste0(i, "\n", title_tmp)) +
          theme_light(base_size = 15) + 
          theme(panel.border = element_rect(fill=NA,color="black", size=1), panel.grid=element_blank())+
          theme(legend.position = "none", plot.title = element_text(hjust = 0.5), aspect.ratio = 1))
  dev.off()
}

for (i in Order){
  markergenes2 = markergenes1[i, ] 
  markergenes2 = as.character(markergenes2)
  Norm_data <- FetchData(sce.SCTnormalize.clean, markergenes2, slot = "data")
  Norm_data$means = rowMeans(Norm_data)
  umap_plot = sce.SCTnormalize.clean@reductions$umap@cell.embeddings %>% as.data.frame()
  new_data = cbind(umap_plot, Norm_data$means)
  pdf(paste0(Figure_path, "/", i, "_marker_average_expression_UMAP_nolabel.pdf"))
  print(ggplot(new_data, aes(x = UMAP_1, y = UMAP_2, colour = Norm_data$means)) + 
          geom_point(size = 0.2, alpha = 1) + 
          scale_colour_gradientn(colours = c("#CCCCCC", "red")) + xlab("") + ylab("") +
          theme_light(base_size = 15) + 
          theme(panel.border = element_rect(fill=NA,color="black", size=1), panel.grid=element_blank())+
          theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
          theme(legend.position = "none", plot.title = element_text(hjust = 0.5), aspect.ratio = 1))
  dev.off()
}

for (i in Order){
  markergenes2 = markergenes1[i, ] 
  markergenes2 = as.character(markergenes2)
  Norm_data <- FetchData(sce.SCTnormalize.clean, markergenes2, slot = "data")
  Norm_data$means = rowMeans(Norm_data)
  tsne_plot = sce.SCTnormalize.clean@reductions$tsne@cell.embeddings %>% as.data.frame()
  new_data = cbind(tsne_plot, Norm_data$means)
  title_tmp = paste0("(", paste0(markergenes2, collapse = ", "),")")
  pdf(paste0(Figure_path, "/", i, "_marker_average_expression_TSNE.pdf"))
  print(ggplot(new_data, aes(x = tSNE_1, y = tSNE_2, colour = Norm_data$means)) + 
          geom_point(size = 0.2, alpha = 1) + 
          scale_colour_gradientn(colours = c("#CCCCCC", "red")) + 
          labs(title = paste0(i, "\n", title_tmp)) +
          theme_light(base_size = 15) + 
          theme(panel.border = element_rect(fill=NA,color="black", size=1), panel.grid=element_blank())+
          theme(legend.position = "none", plot.title = element_text(hjust = 0.5), aspect.ratio = 1))
  dev.off()
}

for (i in Order){
  markergenes2 = markergenes1[i, ] 
  markergenes2 = as.character(markergenes2)
  Norm_data <- FetchData(sce.SCTnormalize.clean, markergenes2, slot = "data")
  Norm_data$means = rowMeans(Norm_data)
  tsne_plot = sce.SCTnormalize.clean@reductions$tsne@cell.embeddings %>% as.data.frame()
  new_data = cbind(tsne_plot, Norm_data$means)
  pdf(paste0(Figure_path, "/", i, "_marker_average_expression_TSNE_nolabel.pdf"))
  print(ggplot(new_data, aes(x = tSNE_1, y = tSNE_2, colour = Norm_data$means)) + 
          geom_point(size = 0.2, alpha = 1) + 
          scale_colour_gradientn(colours = c("#CCCCCC", "red")) + xlab("") + ylab("") +
          theme_light(base_size = 15) + 
          theme(panel.border = element_rect(fill=NA,color="black", size=1), panel.grid=element_blank())+
          theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
          theme(legend.position = "none", plot.title = element_text(hjust = 0.5), aspect.ratio = 1))
  dev.off()
}
