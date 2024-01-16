sample_name = ["T1", "T2", "T3", "T4", "T5", "T6"]
Load_rdata = "figure/sce.SCTnormalize.Rdata"
output_rdata = "doubltremove/sce.SCTnormalize.Doubltremove.Rdata"

import os
if not os.path.exists('doubltremove'):
    os.mkdir('doubltremove')

output_plot = []
for i in sample_name:
    tmp = "doubltremove/DimPlot_"+i+"_doublet_cell.pdf"
    output_plot.append(tmp)

output_file = []
for i in sample_name:
    tmp = "doubltremove/"+i+"_doublet_barcode.txt"
    output_file.append(tmp)

rule all:
    input:
        output_plot,
        output_file,
        "doubltremove/merge_doublet_barcode.txt",
        output_rdata
	
rule DoubltFinder:
    output:
        output_plot,
        output_file
    params:
        sample_name
    script:
        "DoubltFinder_snakemake.R"

rule merge:
    input:
        output_file
    output:
        "doubltremove/merge_doublet_barcode.txt"
    shell:
        "cat {input} > {output}"

rule remove_doublt:
    input:
        "doubltremove/merge_doublet_barcode.txt"
    output:
        output_rdata
    params:
        Load_rdata
    script:
        "remove_doublt_snakemake.R"
