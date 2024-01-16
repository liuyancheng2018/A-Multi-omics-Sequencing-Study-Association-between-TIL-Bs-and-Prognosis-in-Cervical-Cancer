
library(Seurat)
library(SeuratObject)
library(Scissor)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(GSVA)
library(tidyverse)
library(survival)
library(survminer)

load(snakemake@input[[1]])

color_file = read.table(snakemake@input[[2]], sep = "\t", comment.char = "", header = T)
use_color = color_file$color
names(use_color) = color_file$cell_type

Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- "Scissor+"
Scissor_select[infos1$Scissor_neg] <- "Scissor-"
sce.SCTnormalize.clean <- AddMetaData(sce.SCTnormalize.clean, metadata = Scissor_select, col.name = "scissor")
sce.SCTnormalize.clean@meta.data$scissor[sce.SCTnormalize.clean@meta.data$scissor == 0] <- "Background"

pdf(snakemake@output[[1]], width = 8, height = 8)
DimPlot(sce.SCTnormalize.clean, reduction = 'umap', 
        group.by = 'scissor',
        cols = c('grey','royalblue','indianred1'), 
        pt.size = 0.001, order = c("Scissor+","Scissor-"))
dev.off()

pdf(snakemake@output[[2]], width = 8, height = 8)
DimPlot(sce.SCTnormalize.clean, 
        reduction = "umap", 
        group.by = "cell_type", 
        label = T)
dev.off()

pdf(snakemake@output[[3]], width = 8, height = 8)
DimPlot(sce.SCTnormalize.clean, reduction = 'tsne', 
        group.by = 'scissor',
        cols = c('grey','royalblue','indianred1'), 
        pt.size = 0.001, order = c("Scissor+","Scissor-"))
dev.off()

pdf(snakemake@output[[4]], width = 8, height = 8)
DimPlot(sce.SCTnormalize.clean, 
        reduction = "tsne", 
        group.by = "cell_type", 
        label = T)
dev.off()

################################################################################
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

ggsave(snakemake@output[[5]], width = 8, height = 8)
################################################################################
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

ggsave(snakemake@output[[6]], width = 8, height = 8)

################################################################################
Idents(sce.SCTnormalize.clean) = "scissor"
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

ggsave(snakemake@output[[7]], width = 8, height = 8)
################################################################################
Idents(sce.SCTnormalize.clean) = "scissor"
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

ggsave(snakemake@output[[8]], width = 8, height = 8)

################################################################################

data = sce.SCTnormalize.clean@meta.data[, c("cell_type", "scissor")]

proportions <- data %>%
  group_by(cell_type, scissor) %>%
  summarize(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

write.table(proportions, snakemake@output[[9]], sep = "\t", row.names = FALSE)

################################################################################
ggplot(proportions, aes(x = cell_type, y = percentage, fill = scissor)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = use_color) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(title = "", x = "", y = "Percentage") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.ticks.length = unit(.1, "cm"),
    legend.title = element_blank(),
    axis.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1.0)
  )

ggsave(snakemake@output[[10]], width = 6, height = 8)

################################################################################
ggplot(proportions, aes(x = cell_type, y = percentage, fill = scissor)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = use_color) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(title = "", x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.length = unit(.1, "cm"),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  )

ggsave(snakemake@output[[11]], width = 6, height = 8)

#################################################################################
proportions_filtered <- proportions[proportions$scissor == "Scissor+", ]

proportions_filtered$proportion <- round(proportions_filtered$count / sum(proportions_filtered$count) * 100, 2)
proportions_filtered$proportion <- paste(proportions_filtered$proportion, "%", sep = "")
ggplot(proportions_filtered, aes(x = reorder(cell_type, -count), y = count, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste(count, "\n(", proportion, ")")),
            position = position_dodge(width = 0.8), vjust = -0.2, color = "black") +
  scale_fill_manual(values = use_color) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(title = "scissor+", x = "", y = "Number of cells") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.ticks.length = unit(.1, "cm"),
    legend.title = element_blank(),
    axis.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1.0),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(snakemake@output[[12]], width = 6, height = 8)

################################################################################
proportions_filtered <- proportions[proportions$scissor == "Scissor-", ]

proportions_filtered$proportion <- round(proportions_filtered$count / sum(proportions_filtered$count) * 100, 2)
proportions_filtered$proportion <- paste(proportions_filtered$proportion, "%", sep = "")
ggplot(proportions_filtered, aes(x = reorder(cell_type, -count), y = count, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste(count, "\n(", proportion, ")")),
            position = position_dodge(width = 0.8), vjust = -0.2, color = "black") +
  scale_fill_manual(values = use_color) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(title = "scissor-", x = "", y = "Number of cells") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.ticks.length = unit(.1, "cm"),
    legend.title = element_blank(),
    axis.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1.0),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(snakemake@output[[13]], width = 6, height = 8)

################################################################################
proportions_filtered <- proportions[proportions$scissor == "Scissor+", ]

proportions_filtered$proportion <- round(proportions_filtered$count / sum(proportions_filtered$count) * 100, 2)
proportions_filtered$proportion <- paste(proportions_filtered$proportion, "%", sep = "")
ggplot(proportions_filtered, aes(x = reorder(cell_type, -count), y = count, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = use_color) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(title = "", x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.length = unit(.1, "cm"),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  )

ggsave(snakemake@output[[14]], width = 5, height = 8)

################################################################################
proportions_filtered <- proportions[proportions$scissor == "Scissor-", ]

proportions_filtered$proportion <- round(proportions_filtered$count / sum(proportions_filtered$count) * 100, 2)
proportions_filtered$proportion <- paste(proportions_filtered$proportion, "%", sep = "")
ggplot(proportions_filtered, aes(x = reorder(cell_type, -count), y = count, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = use_color) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(title = "", x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.length = unit(.1, "cm"),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  )

ggsave(snakemake@output[[15]], width = 5, height = 8)

identical_1 = "Scissor+"
identical_2 = "Scissor-"

pro =  paste0(identical_1, "_vs_", identical_2)

Idents(sce.SCTnormalize.clean) = "scissor" 

deg = FindMarkers(object = sce.SCTnormalize.clean, 
                  ident.1 = identical_1,
                  ident.2 = identical_2, 
                  test.use='wilcox')

deg$log10pdj = -log10(deg$p_val_adj)

write.table(deg, paste0("result/", pro, "_differential_gene.txt"), sep = "\t", col.names = NA)

dataset <-subset(deg, deg$avg_log2FC!="NA")
dataset <-subset(dataset,dataset$p_val_adj!="NA")

xlim_value = ceiling(max(abs(dataset$avg_log2FC)))
ceiling2 = function(x, base = 50){
  x + ifelse(x%%base == 0, 0, 50-x%%base)
}

ylim_value = ceiling2(max(abs(dataset$log10pdj[!is.infinite(dataset$log10pdj)])))

x = dataset$log10pdj
dataset[, -1] = apply(dataset[, -1], 2, function(x){
  x[is.infinite(x)] <- ylim_value
  x
})

dataset$threshold <- as.factor(ifelse(dataset$p_val_adj < 0.01 & 
                                        abs(dataset$avg_log2FC) >= 1,
                                      ifelse(dataset$avg_log2FC >= 1,'Up','Down'),'Not'))
dataset$gene_name = rownames(dataset)
theme = theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              axis.ticks=element_line(colour="black"),
              plot.margin=unit(c(1,1,1,1),"line"),
              legend.position = 'none',
              legend.key = element_rect(fill = "white"),
              text = element_text(family = "sans",size = 16,face = "bold"))

ggplot(data = dataset, aes(x=avg_log2FC, y=log10pdj, colour=threshold, label = gene_name)) +
  xlab(expression(paste(log[2], "(Fold Change)", sep = "")))+
  ylab(expression(paste(-log[10], "(padj)", sep = ""))) +
  geom_point(size = 2,alpha = 1) + 
  ylim(0, ylim_value) + xlim(-xlim_value, xlim_value) +
  scale_color_manual(values=c("blue", "grey", "red"))+
  geom_vline(xintercept = c(-1, 1), lty = 2,colour="#000000")+
  geom_hline(yintercept = c(2), lty = 2,colour="#000000")+theme+
  geom_text_repel(data = subset(dataset, dataset$p_val_adj < 0.05 & abs(dataset$avg_log2FC) >= 1),
                  aes(x=avg_log2FC, y=log10pdj, label = gene_name), 
                  size = 3, segment.color = "black", show.legend = FALSE)

ggsave(paste0("result/", pro, "_volcanoplot.pdf"), height = 8, width = 6)
write.table(dataset, paste0("result/", pro, "_volcanoplot.txt"), sep = "\t", row.names = F)

survival = phenotype[, c("days_to_last_follow_up", "vital_status")]
rownames(survival) <- gsub("-", ".", rownames(survival))
survival$days_to_last_follow_up = round(survival$days_to_last_follow_up/365, 3)
survival <- rownames_to_column(survival, var = "sample_id")


ident_1_high = rownames(deg[deg$avg_log2FC > 1 & deg$p_val_adj < 0.01,])
ident_2_high = rownames(deg[deg$avg_log2FC < -1 & deg$p_val_adj < 0.01,])

ident_high <- setNames(list(ident_1_high, ident_2_high), c(identical_1, identical_2))

gsva_matrix = gsva(as.matrix(bulk_dataset),
                   ident_high,
                   method = 'ssgsea',
                   kcdf = 'Gaussian',
                   abs.ranking = TRUE)

gsva_result = as.data.frame(t(gsva_matrix))

gsva_result <- rownames_to_column(gsva_result, var = "sample_id")

merge_data = inner_join(survival, gsva_result, by = "sample_id")

for (i in colnames(merge_data)[4:ncol(merge_data)]){
  survival_matrix <- merge_data[, c("sample_id", "vital_status", "days_to_last_follow_up", i)]
  mean_var <- median(survival_matrix[, i])
  survival_matrix$result <- ifelse(survival_matrix[, i] >= mean_var, "high", "low")
  fit <- survfit(Surv(days_to_last_follow_up, vital_status) ~ result, data = survival_matrix)
  pdf(paste0("result/", pro, "_", i, "_gsva_survival_curve.pdf"))
  print(ggsurvplot(fit, data = survival_matrix,
                   linetype = c(1, 1),
                   censor.shape = "",
                   pval = TRUE,
                   conf.int = FALSE,
                   risk.table = FALSE,
                   axes.offset = FALSE,
                   risk.table.col = "strata",
                   surv.median.line = "none",
                   xlab = "Time (years)",
                   ylab = "Overall Survival",
                   title = i,
                   legend = c(0.9, 0.85),
                   legend.title = "",
                   legend.lab = c("high", "low"),
                   ggtheme = theme(text = element_text(size = 12, face = "bold"),
                                   plot.title = element_text(hjust = 0.5),
                                   panel.border = element_blank(),
                                   axis.line = element_line(color="black", size=1),
                                   panel.grid=element_blank(),
                                   panel.background = element_blank(),
                                   legend.key = element_blank(),
                                   aspect.ratio = 1),
                   palette = c("red", "blue")))
  dev.off()
}


for (i in colnames(merge_data)[4:ncol(merge_data)]){
  survival_matrix <- merge_data[, c("sample_id", "vital_status", "days_to_last_follow_up", i)]
  mean_var <- median(survival_matrix[, i])
  survival_matrix$result <- ifelse(survival_matrix[, i] >= mean_var, "high", "low")
  fit <- survfit(Surv(days_to_last_follow_up, vital_status) ~ result, data = survival_matrix)
  pdf(paste0("result/", pro, "_", i, "_gsva_survival_curve_nolabel.pdf"))
  print(ggsurvplot(fit, data = survival_matrix,
                   linetype = c(1, 1),
                   censor.shape = "",
                   pval = FALSE,
                   conf.int = FALSE,
                   risk.table = FALSE,
                   axes.offset = FALSE,
                   risk.table.col = "strata",
                   surv.median.line = "none",
                   xlab = "",
                   ylab = "",
                   title = "",
                   legend = "none",
                   ggtheme = theme(text = element_text(size = 0),
                                   plot.title = element_text(hjust = 0.5),
                                   panel.border = element_blank(),
                                   axis.line = element_line(color="black", size=1),
                                   panel.grid=element_blank(),
                                   panel.background = element_blank(),
                                   legend.key = element_blank(),
                                   aspect.ratio = 1),
                   palette = c("red", "blue")))
  dev.off()
}


################################################################################
theme = theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_line(colour="black"),
              plot.margin=unit(c(1,1,1,1),"line"),
              legend.position = 'none',
              legend.key = element_rect(fill = "white"),
              text = element_text(family = "sans",size = 16,face = "bold"))

ggplot(data = dataset, aes(x=avg_log2FC, y=log10pdj, colour=threshold, label = gene_name)) +
  xlab("")+ylab("") +
  geom_point(size = 2,alpha = 1) + 
  ylim(0, ylim_value) + xlim(-xlim_value, xlim_value) +
  scale_color_manual(values=c("blue", "grey", "red"))+
  geom_vline(xintercept = c(-1, 1), lty = 2,colour="#000000")+
  geom_hline(yintercept = c(2), lty = 2,colour="#000000")+theme

ggsave(paste0("result/", pro, "_volcanoplot_nolabel.pdf"), height = 8, width = 6)
################################################################################
scissor_positive_vs_other <- rep(0, ncol(sce.SCTnormalize.clean))
names(scissor_positive_vs_other) <- colnames(sce.SCTnormalize.clean)
scissor_positive_vs_other <- ifelse(sce.SCTnormalize.clean@meta.data$scissor == "Scissor+", "Scissor+", "Other")
sce.SCTnormalize.clean <- AddMetaData(sce.SCTnormalize.clean, 
                                      metadata = scissor_positive_vs_other, 
                                      col.name = "scissor_positive_vs_other")

identical_1 = "Scissor+"
identical_2 = "Other"

pro =  paste0(identical_1, "_vs_", identical_2)

Idents(sce.SCTnormalize.clean) = "scissor_positive_vs_other" 

deg = FindMarkers(object = sce.SCTnormalize.clean, 
                  ident.1 = identical_1,
                  ident.2 = identical_2, 
                  test.use='wilcox')

deg$log10pdj = -log10(deg$p_val_adj)

write.table(deg, paste0("result/", pro, "_differential_gene.txt"), sep = "\t", col.names = NA)


dataset <-subset(deg, deg$avg_log2FC!="NA")
dataset <-subset(dataset,dataset$p_val_adj!="NA")

xlim_value = ceiling(max(abs(dataset$avg_log2FC)))
ceiling2 = function(x, base = 50){
  x + ifelse(x%%base == 0, 0, 50-x%%base)
}

ylim_value = ceiling2(max(abs(dataset$log10pdj[!is.infinite(dataset$log10pdj)])))

x = dataset$log10pdj
dataset[, -1] = apply(dataset[, -1], 2, function(x){
  x[is.infinite(x)] <- ylim_value
  x
})

dataset$threshold <- as.factor(ifelse(dataset$p_val_adj < 0.01 & 
                                        abs(dataset$avg_log2FC) >= 1,
                                      ifelse(dataset$avg_log2FC >= 1,'Up','Down'),'Not'))
dataset$gene_name = rownames(dataset)
theme = theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              axis.ticks=element_line(colour="black"),
              plot.margin=unit(c(1,1,1,1),"line"),
              legend.position = 'none',
              legend.key = element_rect(fill = "white"),
              text = element_text(family = "sans",size = 16,face = "bold"))

ggplot(data = dataset, aes(x=avg_log2FC, y=log10pdj, colour=threshold, label = gene_name)) +
  xlab(expression(paste(log[2], "(Fold Change)", sep = "")))+
  ylab(expression(paste(-log[10], "(padj)", sep = ""))) +
  geom_point(size = 2,alpha = 1) + 
  ylim(0, ylim_value) + xlim(-xlim_value, xlim_value) +
  scale_color_manual(values=c("blue", "grey", "red"))+
  geom_vline(xintercept = c(-1, 1), lty = 2,colour="#000000")+
  geom_hline(yintercept = c(2), lty = 2,colour="#000000")+theme+
  geom_text_repel(data = subset(dataset, dataset$p_val_adj < 0.05 & abs(dataset$avg_log2FC) >= 1),
                  aes(x=avg_log2FC, y=log10pdj, label = gene_name), 
                  size = 3, segment.color = "black", show.legend = FALSE)

ggsave(paste0("result/", pro, "_volcanoplot.pdf"), height = 8, width = 6)
write.table(dataset, paste0("result/", pro, "_volcanoplot.txt"), sep = "\t", row.names = F)

################################################################################
theme = theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_line(colour="black"),
              plot.margin=unit(c(1,1,1,1),"line"),
              legend.position = 'none',
              legend.key = element_rect(fill = "white"),
              text = element_text(family = "sans",size = 16,face = "bold"))

ggplot(data = dataset, aes(x=avg_log2FC, y=log10pdj, colour=threshold, label = gene_name)) +
  xlab("")+ylab("") +
  geom_point(size = 2,alpha = 1) + 
  ylim(0, ylim_value) + xlim(-xlim_value, xlim_value) +
  scale_color_manual(values=c("blue", "grey", "red"))+
  geom_vline(xintercept = c(-1, 1), lty = 2,colour="#000000")+
  geom_hline(yintercept = c(2), lty = 2,colour="#000000")+theme

ggsave(paste0("result/", pro, "_volcanoplot_nolabel.pdf"), height = 8, width = 6)

ident_1_high = rownames(deg[deg$avg_log2FC > 1 & deg$p_val_adj < 0.01,])
ident_2_high = rownames(deg[deg$avg_log2FC < -1 & deg$p_val_adj < 0.01,])

ident_high <- setNames(list(ident_1_high, ident_2_high), c(identical_1, identical_2))

gsva_matrix = gsva(as.matrix(bulk_dataset),
                   ident_high,
                   method = 'ssgsea',
                   kcdf = 'Gaussian',
                   abs.ranking = TRUE)

gsva_result = as.data.frame(t(gsva_matrix))

gsva_result <- rownames_to_column(gsva_result, var = "sample_id")

merge_data = inner_join(survival, gsva_result, by = "sample_id")

for (i in colnames(merge_data)[4:ncol(merge_data)]){
  survival_matrix <- merge_data[, c("sample_id", "vital_status", "days_to_last_follow_up", i)]
  mean_var <- median(survival_matrix[, i])
  survival_matrix$result <- ifelse(survival_matrix[, i] >= mean_var, "high", "low")
  fit <- survfit(Surv(days_to_last_follow_up, vital_status) ~ result, data = survival_matrix)
  pdf(paste0("result/", pro, "_", i, "_gsva_survival_curve.pdf"))
  print(ggsurvplot(fit, data = survival_matrix,
                   linetype = c(1, 1),
                   censor.shape = "",
                   pval = TRUE,
                   conf.int = FALSE,
                   risk.table = FALSE,
                   axes.offset = FALSE,
                   risk.table.col = "strata",
                   surv.median.line = "none",
                   xlab = "Time (years)",
                   ylab = "Overall Survival",
                   title = i,
                   legend = c(0.9, 0.85),
                   legend.title = "",
                   legend.lab = c("high", "low"),
                   ggtheme = theme(text = element_text(size = 12, face = "bold"),
                                   plot.title = element_text(hjust = 0.5),
                                   panel.border = element_blank(),
                                   axis.line = element_line(color="black", size=1),
                                   panel.grid=element_blank(),
                                   panel.background = element_blank(),
                                   legend.key = element_blank(),
                                   aspect.ratio = 1),
                   palette = c("red", "blue")))
  dev.off()
}


for (i in colnames(merge_data)[4:ncol(merge_data)]){
  survival_matrix <- merge_data[, c("sample_id", "vital_status", "days_to_last_follow_up", i)]
  mean_var <- median(survival_matrix[, i])
  survival_matrix$result <- ifelse(survival_matrix[, i] >= mean_var, "high", "low")
  fit <- survfit(Surv(days_to_last_follow_up, vital_status) ~ result, data = survival_matrix)
  pdf(paste0("result/", pro, "_", i, "_gsva_survival_curve_nolabel.pdf"))
  print(ggsurvplot(fit, data = survival_matrix,
                   linetype = c(1, 1),
                   censor.shape = "",
                   pval = FALSE,
                   conf.int = FALSE,
                   risk.table = FALSE,
                   axes.offset = FALSE,
                   risk.table.col = "strata",
                   surv.median.line = "none",
                   xlab = "",
                   ylab = "",
                   title = "",
                   legend = "none",
                   ggtheme = theme(text = element_text(size = 0),
                                   plot.title = element_text(hjust = 0.5),
                                   panel.border = element_blank(),
                                   axis.line = element_line(color="black", size=1),
                                   panel.grid=element_blank(),
                                   panel.background = element_blank(),
                                   legend.key = element_blank(),
                                   aspect.ratio = 1),
                   palette = c("red", "blue")))
  dev.off()
}

################################################################################
scissor_negative_vs_other <- rep(0, ncol(sce.SCTnormalize.clean))
names(scissor_negative_vs_other) <- colnames(sce.SCTnormalize.clean)
scissor_negative_vs_other <- ifelse(sce.SCTnormalize.clean@meta.data$scissor == "Scissor-", "Scissor-", "Other")
sce.SCTnormalize.clean <- AddMetaData(sce.SCTnormalize.clean, 
                                      metadata = scissor_negative_vs_other, 
                                      col.name = "scissor_negative_vs_other")

identical_1 = "Scissor-"
identical_2 = "Other"

pro =  paste0(identical_1, "_vs_", identical_2)

Idents(sce.SCTnormalize.clean) = "scissor_negative_vs_other" 

deg = FindMarkers(object = sce.SCTnormalize.clean, 
                  ident.1 = identical_1,
                  ident.2 = identical_2, 
                  test.use='wilcox')

deg$log10pdj = -log10(deg$p_val_adj)

write.table(deg, paste0("result/", pro, "_differential_gene.txt"), sep = "\t", col.names = NA)


dataset <-subset(deg, deg$avg_log2FC!="NA")
dataset <-subset(dataset,dataset$p_val_adj!="NA")

xlim_value = ceiling(max(abs(dataset$avg_log2FC)))
ceiling2 = function(x, base = 50){
  x + ifelse(x%%base == 0, 0, 50-x%%base)
}

ylim_value = ceiling2(max(abs(dataset$log10pdj[!is.infinite(dataset$log10pdj)])))

x = dataset$log10pdj
dataset[, -1] = apply(dataset[, -1], 2, function(x){
  x[is.infinite(x)] <- ylim_value
  x
})

dataset$threshold <- as.factor(ifelse(dataset$p_val_adj < 0.01 & 
                                        abs(dataset$avg_log2FC) >= 1,
                                      ifelse(dataset$avg_log2FC >= 1,'Up','Down'),'Not'))
dataset$gene_name = rownames(dataset)
theme = theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              axis.ticks=element_line(colour="black"),
              plot.margin=unit(c(1,1,1,1),"line"),
              legend.position = 'none',
              legend.key = element_rect(fill = "white"),
              text = element_text(family = "sans",size = 16,face = "bold"))

ggplot(data = dataset, aes(x=avg_log2FC, y=log10pdj, colour=threshold, label = gene_name)) +
  xlab(expression(paste(log[2], "(Fold Change)", sep = "")))+
  ylab(expression(paste(-log[10], "(padj)", sep = ""))) +
  geom_point(size = 2,alpha = 1) + 
  ylim(0, ylim_value) + xlim(-xlim_value, xlim_value) +
  scale_color_manual(values=c("blue", "grey", "red"))+
  geom_vline(xintercept = c(-1, 1), lty = 2,colour="#000000")+
  geom_hline(yintercept = c(2), lty = 2,colour="#000000")+theme+
  geom_text_repel(data = subset(dataset, dataset$p_val_adj < 0.05 & abs(dataset$avg_log2FC) >= 1),
                  aes(x=avg_log2FC, y=log10pdj, label = gene_name), 
                  size = 3, segment.color = "black", show.legend = FALSE)

ggsave(paste0("result/", pro, "_volcanoplot.pdf"), height = 8, width = 6)
write.table(dataset, paste0("result/", pro, "_volcanoplot.txt"), sep = "\t", row.names = F)

################################################################################
theme = theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_line(colour="black"),
              plot.margin=unit(c(1,1,1,1),"line"),
              legend.position = 'none',
              legend.key = element_rect(fill = "white"),
              text = element_text(family = "sans",size = 16,face = "bold"))

ggplot(data = dataset, aes(x=avg_log2FC, y=log10pdj, colour=threshold, label = gene_name)) +
  xlab("")+ylab("") +
  geom_point(size = 2,alpha = 1) + 
  ylim(0, ylim_value) + xlim(-xlim_value, xlim_value) +
  scale_color_manual(values=c("blue", "grey", "red"))+
  geom_vline(xintercept = c(-1, 1), lty = 2,colour="#000000")+
  geom_hline(yintercept = c(2), lty = 2,colour="#000000")+theme

ggsave(paste0("result/", pro, "_volcanoplot_nolabel.pdf"), height = 8, width = 6)

ident_1_high = rownames(deg[deg$avg_log2FC > 1 & deg$p_val_adj < 0.01,])
ident_2_high = rownames(deg[deg$avg_log2FC < -1 & deg$p_val_adj < 0.01,])

ident_high <- setNames(list(ident_1_high, ident_2_high), c(identical_1, identical_2))

gsva_matrix = gsva(as.matrix(bulk_dataset),
                   ident_high,
                   method = 'ssgsea',
                   kcdf = 'Gaussian',
                   abs.ranking = TRUE)

gsva_result = as.data.frame(t(gsva_matrix))

gsva_result <- rownames_to_column(gsva_result, var = "sample_id")

merge_data = inner_join(survival, gsva_result, by = "sample_id")

for (i in colnames(merge_data)[4:ncol(merge_data)]){
  survival_matrix <- merge_data[, c("sample_id", "vital_status", "days_to_last_follow_up", i)]
  mean_var <- median(survival_matrix[, i])
  survival_matrix$result <- ifelse(survival_matrix[, i] >= mean_var, "high", "low")
  fit <- survfit(Surv(days_to_last_follow_up, vital_status) ~ result, data = survival_matrix)
  pdf(paste0("result/", pro, "_", i, "_gsva_survival_curve.pdf"))
  print(ggsurvplot(fit, data = survival_matrix,
                   linetype = c(1, 1),
                   censor.shape = "",
                   pval = TRUE,
                   conf.int = FALSE,
                   risk.table = FALSE,
                   axes.offset = FALSE,
                   risk.table.col = "strata",
                   surv.median.line = "none",
                   xlab = "Time (years)",
                   ylab = "Overall Survival",
                   title = i,
                   legend = c(0.9, 0.85),
                   legend.title = "",
                   legend.lab = c("high", "low"),
                   ggtheme = theme(text = element_text(size = 12, face = "bold"),
                                   plot.title = element_text(hjust = 0.5),
                                   panel.border = element_blank(),
                                   axis.line = element_line(color="black", size=1),
                                   panel.grid=element_blank(),
                                   panel.background = element_blank(),
                                   legend.key = element_blank(),
                                   aspect.ratio = 1),
                   palette = c("red", "blue")))
  dev.off()
}


for (i in colnames(merge_data)[4:ncol(merge_data)]){
  survival_matrix <- merge_data[, c("sample_id", "vital_status", "days_to_last_follow_up", i)]
  mean_var <- median(survival_matrix[, i])
  survival_matrix$result <- ifelse(survival_matrix[, i] >= mean_var, "high", "low")
  fit <- survfit(Surv(days_to_last_follow_up, vital_status) ~ result, data = survival_matrix)
  pdf(paste0("result/", pro, "_", i, "_gsva_survival_curve_nolabel.pdf"))
  print(ggsurvplot(fit, data = survival_matrix,
                   linetype = c(1, 1),
                   censor.shape = "",
                   pval = FALSE,
                   conf.int = FALSE,
                   risk.table = FALSE,
                   axes.offset = FALSE,
                   risk.table.col = "strata",
                   surv.median.line = "none",
                   xlab = "",
                   ylab = "",
                   title = "",
                   legend = "none",
                   ggtheme = theme(text = element_text(size = 0),
                                   plot.title = element_text(hjust = 0.5),
                                   panel.border = element_blank(),
                                   axis.line = element_line(color="black", size=1),
                                   panel.grid=element_blank(),
                                   panel.background = element_blank(),
                                   legend.key = element_blank(),
                                   aspect.ratio = 1),
                   palette = c("red", "blue")))
  dev.off()
}

proportions_filtered <- proportions[proportions$scissor == "Scissor+", ]

ggplot(proportions_filtered, aes(x = reorder(cell_type, -percentage), y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(round(percentage, 2), "%")),
            position = position_dodge(width = 0.8), vjust = -0.2, color = "black") +
  scale_fill_manual(values = use_color) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(title = "", x = "", y = "Scossor+ cell number / Total cell number") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.ticks.length = unit(.1, "cm"),
    legend.title = element_blank(),
    axis.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1.0),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(snakemake@output[[16]], width = 6, height = 8)

proportions_filtered <- proportions[proportions$scissor == "Scissor+", ]

ggplot(proportions_filtered, aes(x = reorder(cell_type, -percentage), y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = use_color) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(title = "", x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.length = unit(.1, "cm"),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  )

ggsave(snakemake@output[[17]], width = 5, height = 8)


proportions_filtered <- proportions[proportions$scissor == "Scissor-", ]

ggplot(proportions_filtered, aes(x = reorder(cell_type, -percentage), y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(round(percentage, 2), "%")),
            position = position_dodge(width = 0.8), vjust = -0.2, color = "black") +
  scale_fill_manual(values = use_color) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(title = "", x = "", y = "Scossor- cell number / Total cell number") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.ticks.length = unit(.1, "cm"),
    legend.title = element_blank(),
    axis.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1.0),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(snakemake@output[[18]], width = 6, height = 8)

proportions_filtered <- proportions[proportions$scissor == "Scissor-", ]

ggplot(proportions_filtered, aes(x = reorder(cell_type, -percentage), y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = use_color) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(title = "", x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.length = unit(.1, "cm"),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none"
  )

ggsave(snakemake@output[[19]], width = 5, height = 8)




