pFilter = 0.1

rule all:
    input:
        "result/multi_Cox.txt",
        "result/multi_SigExp.txt",
        "result/multi_Cox_forest.pdf"

rule multicox:
    input:
        "result/lasso_SigExp.txt"
    output:
        "result/multi_Cox.txt",
        "result/multi_SigExp.txt"
    params:
        pFilter
    script:
        "multicox.R"
rule forest:
    input:
        "result/multi_Cox.txt"
    output:
        "result/multi_Cox_forest.pdf"
    script:
        "forest.R"
