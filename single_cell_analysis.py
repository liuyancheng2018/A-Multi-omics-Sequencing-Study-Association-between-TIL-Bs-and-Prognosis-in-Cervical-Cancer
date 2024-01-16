Dimension_value = 30
Resolution_value = 0.6
Seu_obj = "sce.SCTnormalize.clean"
Cell_type = "all_cells"

plot = 'plot_Dimension_value_'+str(Dimension_value)+'_Resolution_value_'+str(Resolution_value)

import pandas as pd

cell_type = pd.read_csv('all_cell_type.txt', delimiter='\t')
Sample = cell_type['cell_type'].tolist()

cluster_cell_type = pd.read_csv('all_cluster_cell_type.txt', delimiter='\t')
cluster2celltype = cluster_cell_type['cell_type'].tolist()

Load_Rdata = "heatmap_dimention_value_30_resolution_value_0.6/sce.SCTnormalize.Doubltremove.all.30.0.6.Rdata"
output_file_1 = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_celltype_umap.pdf"
output_file_2 = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_celltype_tsne.pdf"
output_file_5 = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_sample_umap.pdf"
output_file_6 = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_sample_tsne.pdf"
output_file_7 = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_celltype_tsne_good.pdf"
output_file_8 = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_celltype_tsne_good_nolabel.pdf"
output_file_9 = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_celltype_umap_good.pdf"
output_file_10 = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_celltype_umap_good_nolabel.pdf"
output_file_24 = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_celltype_umap_good_black.pdf"
output_file_28 = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_celltype_tsne_good_black.pdf"
output_rdata = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_celltype.Rdata"
barplot_file = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_cell_ratio_bar.pdf"
barplot_file_nolabel = plot+"/"+Cell_type+"_dimension_"+str(Dimension_value)+"_resolution_"+str(Resolution_value)+"_cell_ratio_bar_nolabel.pdf"
Order = ["T1", "T2", "T3", "T4", "T5", "T6"]


cellnumber_statistics_files = []
for i in Order:
    tmp = plot+"/"+i+"_cell_number_statistics.txt"
    cellnumber_statistics_files.append(tmp)

sample_cell_number_statistics = plot+"/sample_cell_number_statistics.txt"
celltype_cell_number_statistics = plot+"/celltype_cell_number_statistics.txt"
merge_statistics = plot+"/merge_statistics.txt"

celltype_compare_p_value = plot+"/celltype_compare_p_value.txt",
cell_ratio_compare_boxplot = plot+"/cell_ratio_compare_boxplot.pdf",
cell_ratio_compare_boxplot_nolabel = plot+"/cell_ratio_compare_boxplot_nolabel.pdf",
cell_ratio_compare_barplot = plot+"/cell_ratio_compare_barplot.pdf",
cell_ratio_compare_barplot_nolabel = plot+"/cell_ratio_compare_barplot_nolabel.pdf"

color_file = pd.read_csv('color.txt', delimiter='\t')
color_file_cell_type = color_file['cell_type'].tolist()
lagend_color = []
for i in color_file_cell_type:
    tmp = plot+"/label_"+i+".pdf"
    lagend_color.append(tmp)

color_file_2 = pd.read_csv('sample_color.txt', delimiter='\t')
color_file_cell_type_2 = color_file_2['sample'].tolist()
lagend_color_2 = []
for i in color_file_cell_type_2:
    tmp = plot+"/label_"+i+".pdf"
    lagend_color_2.append(tmp)

rule all:
    input:
        output_file_1,
        output_file_2,
        output_file_5,
        output_file_6,
        output_file_7,
        output_file_8,
        output_file_9,
        output_file_10,
        output_rdata,
        output_file_24,
        output_file_28,
        lagend_color,
        barplot_file,
        barplot_file_nolabel,
        sample_cell_number_statistics,
        celltype_cell_number_statistics,
        cellnumber_statistics_files,
        merge_statistics,
        celltype_compare_p_value,
        cell_ratio_compare_boxplot,
        cell_ratio_compare_boxplot_nolabel,
        cell_ratio_compare_barplot,
        cell_ratio_compare_barplot_nolabel,
        lagend_color_2

rule prepare_4:
    input:
        "color.txt"
    params:
        Load_Rdata,
        cluster2celltype,
        plot
    output:
        output_file_1,
        output_file_2,
        output_file_5,
        output_file_6,
        output_file_7,
        output_file_8,
        output_file_9,
        output_file_10,
        output_rdata,
        output_file_24,
        output_file_28,
        lagend_color
    script:
        "all_single_prepare_4.R"

rule barplot:
    input:
        "color.txt",
        output_rdata
    output:
        barplot_file,
        barplot_file_nolabel,
        sample_cell_number_statistics,
        celltype_cell_number_statistics,
        cellnumber_statistics_files
    params:
        Seu_obj,
        Order,
        plot
    script:
        "barplot_snakemake.R"

rule merge_statistics:
    input:
        cellnumber_statistics_files
    output:
        merge_statistics
    shell:
        "cat {input} > {output}"

rule cell_ratio_compare:
    input:
        merge_statistics,
        "sample_color.txt"
    output:
        celltype_compare_p_value,
        cell_ratio_compare_boxplot,
        cell_ratio_compare_boxplot_nolabel,
        cell_ratio_compare_barplot,
        cell_ratio_compare_barplot_nolabel,
        lagend_color_2
    params:
        Sample,
        plot
    script:
        "cell_ratio_compare_boxplot_snakemake.R"
