
rule all:
    input:
        "result/cluster_group2.txt"

rule consensuscluster:
    input:
        "uniSigExp.txt"
    output:
        "result/cluster_group2.txt"
    script:
        "consensuscluster.R"
