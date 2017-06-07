library('MASS')

args <- commandArgs(trailingOnly = TRUE)

top <- as.integer(args[1]) #Top X percent of calls

data <- read.table(args[2], header=FALSE, sep = ";")

data <- data[data$V15<0,]
data.sort <- data[order(data$V15),]

#Create observed and background sets
index <- as.integer(nrow(data.sort)*top/100)
data.obs <- data.sort[1:index,] 
data.back <- data.sort[(index+1):nrow(data.sort),]

energies.obs <- data.obs[,15]
energies.back <- data.back[,15]

energies.obs <- energies.obs*(-1)
energies.back <- energies.back*(-1)

#Estimate parameters of background gamma distribution
op <- options(digits = 3)
set.seed(123)
params <- fitdistr(energies.back, "gamma", lower = 0.001)

sh <- params$estimate["shape"]
rt <- params$estimate["rate"]

#Compute p-values and FDR
intarna.pvalues <- 1-pgamma(energies.obs, shape=sh, rate=rt)
intarna.pvalues.corr <- p.adjust(intarna.pvalues, method = "BH")
data.obs[,"pvalues"] <- intarna.pvalues
data.obs[,"FDR"] <- intarna.pvalues.corr

write.table(file=paste(args[2],".pvalues", sep = ""), data.obs, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#warnings()