
rule all:
    input:
        "result/uniForest.pdf",

rule forest_plot:
    input:
        "result/uniCox.txt",
    output:
        "result/uniForest.pdf",
    script:
        "uni_forest.R"
