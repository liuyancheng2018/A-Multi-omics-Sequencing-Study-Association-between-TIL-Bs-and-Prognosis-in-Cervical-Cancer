library(survival)
library(survminer)
library(glmnet)

pFilter = snakemake@params[[1]]

rt = read.table(snakemake@input[[1]], header = T, sep = "\t", check.names = F, row.names = 1)
multiCox = coxph(Surv(days_to_last_follow_up, vital_status) ~ ., data = rt)
# multiCox = step(multiCox, direction = "bath") 
multiCoxSum = summary(multiCox)
outTab = data.frame()
outTab = cbind(
  coef = multiCoxSum$coefficients[, "coef"],
  HR = multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L = multiCoxSum$conf.int[, "lower .95"],
  HR.95H = multiCoxSum$conf.int[, "upper .95"],
  pvalue = multiCoxSum$coefficients[, "Pr(>|z|)"]
)
outTab = cbind(id = row.names(outTab), outTab)
outTab = as.data.frame(gsub("'", "", outTab))
outTab <- outTab[outTab$pvalue <= pFilter, ]
write.table(outTab, file = snakemake@output[[1]], sep = "\t", row.names = F, quote = F)

multicoxGene = c("days_to_last_follow_up", "vital_status", rownames(outTab))
multicoxExp = rt[, multicoxGene]
multicoxExp = cbind(id = row.names(multicoxExp), multicoxExp)
write.table(multicoxExp, file = snakemake@output[[2]], sep = "\t", row.names = F, quote = F)
