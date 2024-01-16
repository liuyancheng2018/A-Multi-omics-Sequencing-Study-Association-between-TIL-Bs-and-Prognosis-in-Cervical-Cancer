library(ConsensusClusterPlus)

data = read.table(snakemake@input[[1]], sep = "\t", header = T, row.names = 1)
data <- t(data[, -(1:2)])

maxK = 9
result = ConsensusClusterPlus(data,
                              maxK = maxK,
                              reps = 50,
                              pItem = 0.8,
                              pFeature = 1,
                              title = "result/ConsensusCluster_pit",
                              clusterAlg = "km",
                              distance = "euclidean",
                              seed = 2013,
                              plot = "png")

cluster = as.data.frame(result[[2]][3])
write.table(cluster, file = snakemake@output[[1]], sep = "\t", quote = F)