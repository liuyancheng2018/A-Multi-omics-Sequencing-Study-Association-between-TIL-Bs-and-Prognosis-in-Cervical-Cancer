
single_cell_rdata = "all_cells_dimension_30_resolution_0.6_celltype.Rdata"

rule all:
    input:
        "result/Scissor_survival.RData",
        "result/Scissor_result.RData",
        "result/DimPlot_scissor_umap.pdf",
        "result/DimPlot_cell_type_umap.pdf",
        "result/DimPlot_scissor_tsne.pdf",
        "result/DimPlot_cell_type_tsne.pdf",
        "result/DimPlot_cell_type_umap_nolabel.pdf",
        "result/DimPlot_cell_type_tsne_nolabel.pdf",
        "result/DimPlot_scissor_umap_nolabel.pdf",
        "result/DimPlot_scissor_tsne_nolabel.pdf",
        "result/scissor_proportions.txt",
        "result/scissor_proportions_barplot.pdf",
        "result/scissor_proportions_barplot_nolabel.pdf",
        "result/scissor_positive_proportions_barplot.pdf",
        "result/scissor_negative_proportions_barplot.pdf",
        "result/scissor_positive_proportions_barplot_nolabel.pdf",
        "result/scissor_negative_proportions_barplot_nolabel.pdf",
        "result/scissor_positive_vs_total_number_barplot.pdf",
        "result/scissor_positive_vs_total_number_barplot_nolabel.pdf",
        "result/scissor_negative_vs_total_number_barplot.pdf",
        "result/scissor_negative_vs_total_number_barplot_nolabel.pdf"

rule scissor:
    input:
        single_cell_rdata,
        "merge_fpkm_final_tumoronly_clinical_have.tsv",
        "clinical_matrix_have.tsv"
    output:
        "result/Scissor_survival.RData",
        "result/Scissor_result.RData"
    script:
        "Scissor.R"

rule scissor_plot:
    input:
        "result/Scissor_result.RData",
        "color.txt"
    output:
        "result/DimPlot_scissor_umap.pdf",
        "result/DimPlot_cell_type_umap.pdf",
        "result/DimPlot_scissor_tsne.pdf",
        "result/DimPlot_cell_type_tsne.pdf",
        "result/DimPlot_cell_type_umap_nolabel.pdf",
        "result/DimPlot_cell_type_tsne_nolabel.pdf",
        "result/DimPlot_scissor_umap_nolabel.pdf",
        "result/DimPlot_scissor_tsne_nolabel.pdf",
        "result/scissor_proportions.txt",
        "result/scissor_proportions_barplot.pdf",
        "result/scissor_proportions_barplot_nolabel.pdf",
        "result/scissor_positive_proportions_barplot.pdf",
        "result/scissor_negative_proportions_barplot.pdf",
        "result/scissor_positive_proportions_barplot_nolabel.pdf",
        "result/scissor_negative_proportions_barplot_nolabel.pdf",
        "result/scissor_positive_vs_total_number_barplot.pdf",
        "result/scissor_positive_vs_total_number_barplot_nolabel.pdf",
        "result/scissor_negative_vs_total_number_barplot.pdf",
        "result/scissor_negative_vs_total_number_barplot_nolabel.pdf"
        
    script:
        "Scissor_plot.R"
