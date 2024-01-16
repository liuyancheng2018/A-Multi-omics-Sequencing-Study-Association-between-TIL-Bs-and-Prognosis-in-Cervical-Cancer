library(dplyr)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(circlize)

cluster = read.table(snakemake@input[[1]], header = T, sep = "\t")
cluster$sample_id <- rownames(cluster)
clinical = read.table(snakemake@input[[2]], sep = "\t", header = T)
colnames(clinical)[1] = "sample_id"
clinical$sample_id <- gsub("-", ".", clinical$sample_id)
new_data = inner_join(clinical, cluster, by = "sample_id")
new_data$days_to_last_follow_up = new_data$days_to_last_follow_up/365
diff = survdiff(Surv(days_to_last_follow_up, vital_status) ~ consensusClass, data = new_data)
pValue = 1-pchisq(diff$chisq, df = 1)
pValue = signif(pValue, 4)
pValue = format(pValue, scientific = TRUE)
fit = survfit(Surv(days_to_last_follow_up, vital_status) ~ consensusClass, data = new_data)
pdf(file = snakemake@output[[1]], width = 6, height = 6)
ggsurvplot(fit, data = new_data,
           linetype = c(1, 1),
           censor.shape = "",
           pval = TRUE,
           conf.int = FALSE,
           risk.table = FALSE,
           axes.offset = FALSE,
           risk.table.col = "strata",
           surv.median.line = "none",
           xlab = "Time (year)",
           ylab = "Overall Survival rate",
           title = paste("Survival curve (p = ", pValue, ")", sep = ""),
           legend = c(0.9, 0.85),
           legend.title = "",
           legend.lab = c("Cluster1", "Cluster2"),
           ggtheme = theme(text = element_text(size = 18, face = "bold"),
                           plot.title = element_text(hjust = 0.5),
                           panel.border = element_blank(),
                           axis.line = element_line(color="black", size=1),
                           panel.grid=element_blank(),
                           panel.background = element_blank(),
                           legend.key = element_blank(),
                           aspect.ratio = 1),
           palette = c("blue", "red"))
dev.off()

pdf(file = snakemake@output[[2]], width = 6, height = 6)
ggsurvplot(fit, data = new_data,
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
           ggtheme = theme(text = element_text(size = 0, face = "bold"),
                           plot.title = element_text(hjust = 0.5),
                           panel.border = element_blank(),
                           axis.line = element_line(color="black", size=1),
                           panel.grid=element_blank(),
                           panel.background = element_blank(),
                           legend.key = element_blank(),
                           aspect.ratio = 1),
           palette = c("blue", "red"))
dev.off()

##########################################################################################
matrix = read.table(snakemake@input[[3]], sep = "\t", header = T)
matrix <- matrix[, -c(2, 3)]
colnames(matrix)[1] = "sample_id"

new_data_2 = inner_join(new_data, matrix, by = "sample_id")
new_data_2 <- new_data_2[order(new_data_2[, "consensusClass"]), ]
new_data_2$age_at_index <- ifelse(new_data_2$age_at_index >= 50, ">=50", "<50")
new_data_3 <- new_data_2[, -(1:10)]
new_data_3 = t(new_data_3)
new_data_4 <- t(apply(new_data_3, 1, scale))

#########################################################################################
selected_columns <- new_data_2[, c(2, 5:9)]
colnames(selected_columns) = c("Age", "M", "N", "T", "Stage", "Pathology")
selected_columns_2 <- new_data_2[, c(2, 5:10)]
colnames(selected_columns_2) = c("Age", "M", "N", "T", "Stage", "Pathology", "Cluster")
Order = colnames(selected_columns)

p.value.vector  = c()
for (i in Order){
  data = selected_columns_2[c("Cluster", i)]
  colnames(data) = c("Cluster", "clinical")
  data = data[(data[, "clinical"] != "Unknown"),]
  tableStat = table(data)
  stat = chisq.test(tableStat)
  p.value.vector = append(p.value.vector, stat$p.value)
  }
  
p.value.df = data.frame(Cluster = Order, p.value = p.value.vector)
p.value.df$significance = ifelse(p.value.df$p.value < 0.0001, "****",
                                 ifelse(p.value.df$p.value < 0.001, "***",
                                        ifelse(p.value.df$p.value < 0.01, "**",
                                               ifelse(p.value.df$p.value < 0.05, "*", ""))))

p.value.df$label = paste0("p=", round(p.value.df$p.value, 3))
p.value.df$label = ifelse(p.value.df$p.value < 0.05, p.value.df$label, "")

write.table(p.value.df, snakemake@output[[3]], sep = "\t", row.names = F)


col_fun = colorRamp2(c(-4, 0, 4), c("#00CCFF", "white", "#FF0000"))

ha = HeatmapAnnotation(
  T = as.factor(new_data_2$ajcc_pathologic_t),
  N = as.factor(new_data_2$ajcc_pathologic_n),
  M = as.factor(new_data_2$ajcc_pathologic_m),
  Age = as.factor(new_data_2$age_at_index),
  Pathology = as.factor(new_data_2$primary_diagnosis),
  Stage = as.factor(new_data_2$figo_stage),
  Cluster = as.factor(new_data_2$consensusClass),
  border = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 12),
  col = list(Cluster = c("1" = "#2ECC46", "2" = "#DA61D8"),
             Stage = c("I" = "#F62B93",
                       "II" = "#EA3F1D",
                       "III" = "#9CC816",
                       "IV" = "#CF6D3A",
                       "Unknown" = "#6D7D74"),
             Pathology = c("Squamous cell carcinoma" = "#FFFF00",
                           "Adenosquamous carcinoma" = "#8B6508",
                           "Adenocarcinoma" = "#6495ED"),
             T = c("T1" = "#C1EFFB",
                   "T2" = "#977EE5",
                   "T3" = "#F890AC",
                   "T4" = "#96F5C8",
                   "Tis" = "#BBD5FC",
                   "TX" = "#F382C1",
                   "Unknown" = "#6D7D74"),
             N = c("N0" = "#89FB86",
                   "N1" = "#89CDEC",
                   "NX" = "#F6A896",
                   "Unknown" = "#6D7D74"),
             M = c("M0" = "#A05D74",
                   "M1" = "#B2F3D9",
                   "MX" = "#8A2BE2",
                   "Unknown" = "#6D7D74"),
             Age = c("<50" = "#C5EAF8",
                     ">=50" = "#BB6641")))
pdf(snakemake@output[[4]], width = 8, height = 5)
Heatmap(new_data_4, name = "Expression",
        column_title = NULL,
        show_row_names = TRUE,
        cluster_columns = FALSE,
        show_column_dend = FALSE, 
        show_row_dend = FALSE, 
        col = col_fun,
        top_annotation = ha,
        column_split = new_data_2$consensusClass)
dev.off()


ha = HeatmapAnnotation(show_legend = FALSE,
                       show_annotation_name = FALSE,
                       T = as.factor(new_data_2$ajcc_pathologic_t),
                       N = as.factor(new_data_2$ajcc_pathologic_n),
                       M = as.factor(new_data_2$ajcc_pathologic_m),
                       Age = as.factor(new_data_2$age_at_index),
                       Pathology = as.factor(new_data_2$primary_diagnosis),
                       Stage = as.factor(new_data_2$figo_stage),
                       Cluster = as.factor(new_data_2$consensusClass),
                       border = TRUE,
                       annotation_name_side = "left",
                       annotation_name_gp = gpar(fontsize = 12),
                       col = list(Cluster = c("1" = "#2ECC46", "2" = "#DA61D8"),
                                  Stage = c("I" = "#F62B93",
                                            "II" = "#EA3F1D",
                                            "III" = "#9CC816",
                                            "IV" = "#CF6D3A",
                                            "Unknown" = "#6D7D74"),
                                  Pathology = c("Squamous cell carcinoma" = "#FFFF00",
                                                "Adenosquamous carcinoma" = "#8B6508",
                                                "Adenocarcinoma" = "#6495ED"),
                                  T = c("T1" = "#C1EFFB",
                                        "T2" = "#977EE5",
                                        "T3" = "#F890AC",
                                        "T4" = "#96F5C8",
                                        "Tis" = "#BBD5FC",
                                        "TX" = "#F382C1",
                                        "Unknown" = "#6D7D74"),
                                  N = c("N0" = "#89FB86",
                                        "N1" = "#89CDEC",
                                        "NX" = "#F6A896",
                                        "Unknown" = "#6D7D74"),
                                  M = c("M0" = "#A05D74",
                                        "M1" = "#B2F3D9",
                                        "MX" = "#8A2BE2",
                                        "Unknown" = "#6D7D74"),
                                  Age = c("<50" = "#C5EAF8",
                                          ">=50" = "#BB6641")))

pdf(snakemake@output[[5]], width = 5, height = 8)
Heatmap(new_data_4, name = "Expression",
        show_heatmap_legend = FALSE,
        column_title = NULL,
        show_row_names = FALSE,
        cluster_columns = FALSE,
        show_column_dend = FALSE, 
        show_row_dend = FALSE, 
        col = col_fun,
        top_annotation = ha,
        column_split = new_data_2$consensusClass)
dev.off()
