Load_rdata = "plot_Dimension_value_30_Resolution_value_0.6/all_cells_dimension_30_resolution_0.6_celltype.Rdata"
Order = ["Epithelial_cells", "Endothelial_cells", "Fibroblasts", "Monocyte_Macrophage", "Neutrophils", "Dendritic_cells", "Mast_cells", "NK_T_cells", "B_cells", "Plasma_cells"]
Seu_obj = "sce.SCTnormalize.clean"
Dimension_value = 30
Resolution_value = 0.6
Cell_type = "All_cells"

plot = 'plot_Dimension_value_'+str(Dimension_value)+'_Resolution_value_'+str(Resolution_value)

marker_average_expression_UMAP = []
for i in Order:
    names = names = plot+"/"+i+"_marker_average_expression_UMAP.pdf"
    marker_average_expression_UMAP.append(names)

marker_average_expression_UMAP_nolabel = []
for i in Order:
    names = names = plot+"/"+i+"_marker_average_expression_UMAP_nolabel.pdf"
    marker_average_expression_UMAP_nolabel.append(names)
    
marker_average_expression_TSNE = []
for i in Order:
    names = names = plot+"/"+i+"_marker_average_expression_TSNE.pdf"
    marker_average_expression_TSNE.append(names)

marker_average_expression_TSNE_nolabel = []
for i in Order:
    names = names = plot+"/"+i+"_marker_average_expression_TSNE_nolabel.pdf"
    marker_average_expression_TSNE_nolabel.append(names)
    
Bubble_Plot = plot+"/Bubble_Plot.pdf"
Bubble_plot_nolabel = plot+"/Bubble_plot_nolabel.pdf"
Bubble_Plot_rev = plot+"/Bubble_Plot_rev.pdf"
Bubble_plot_rev_nolabel = plot+"/Bubble_plot_rev_nolabel.pdf"
Vln_Plot = plot+"/Vln_Plot.pdf"
Vln_plot_nolabel = plot+"/Vln_plot_nolabel.pdf"
Vln_plot_rev_1 = plot+"/Vln_plot_rev_1.pdf"
Vln_plot_rev_1_nolabel = plot+"/Vln_plot_rev_1_nolabel.pdf"
Vln_plot_rev_2 = plot+"/Vln_plot_rev_2.pdf"
Vln_plot_rev_2_nolabel = plot+"/Vln_plot_rev_2_nolabel.pdf"

rule all:
    input:
        Bubble_Plot,
        Bubble_plot_nolabel,
        Bubble_Plot_rev,
        Bubble_plot_rev_nolabel,
        Vln_Plot,
        Vln_plot_nolabel,
        Vln_plot_rev_1,
        Vln_plot_rev_1_nolabel,
        Vln_plot_rev_2,
        Vln_plot_rev_2_nolabel,
        marker_average_expression_UMAP,
        marker_average_expression_UMAP_nolabel,
        marker_average_expression_TSNE,
        marker_average_expression_TSNE_nolabel
	
rule bubbleplot:
    input:
        "all_markergenes.txt"
    output:
        Bubble_Plot,
        Bubble_plot_nolabel,
        Bubble_Plot_rev,
        Bubble_plot_rev_nolabel,
        Vln_Plot,
        Vln_plot_nolabel,
        Vln_plot_rev_1,
        Vln_plot_rev_1_nolabel,
        Vln_plot_rev_2,
        Vln_plot_rev_2_nolabel,
        marker_average_expression_UMAP,
        marker_average_expression_UMAP_nolabel,
        marker_average_expression_TSNE,
        marker_average_expression_TSNE_nolabel        
    params:
        Load_rdata,
        Seu_obj,
        Order,
        plot
    script:
        "all_marker_perform.R"
