
rule all:
    input:
        "result/cluster_group2_survival.pdf",
        "result/cluster_group2_survival_nolabel.pdf",
        "result/cluster_group2_clinical_pvalue.txt",
        "result/cluster_group2_heatmap.pdf",
        "result/cluster_group2_heatmap_nolabel.pdf"

rule consensuscluster_plot:
    input:
        "result/cluster_group2.txt",
        "clinical_matrix_have.tsv",
        "uniSigExp.txt"
    output:
        "result/cluster_group2_survival.pdf",
        "result/cluster_group2_survival_nolabel.pdf",
        "result/cluster_group2_clinical_pvalue.txt",
        "result/cluster_group2_heatmap.pdf",
        "result/cluster_group2_heatmap_nolabel.pdf"
    script:
        "consensuscluster_plot.R"
