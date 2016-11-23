library(evd)

args <- commandArgs(trailingOnly = TRUE)

#data <- as.matrix(read.table(args[1], header=FALSE, sep = "\t"))
data <- read.table(args[1], header=FALSE, sep = "\t")

data.obs <- data[grep("Random", data[,1], invert = TRUE), ] #there might be a faster way
data.back <- data[grep("Random", data[,1]), ]

data.obs.sort <- data.obs[order(data.obs$V7),]
data.back.sort <- data.back[order(data.back$V7),]

energies.obs <- data.obs.sort[,7]
energies.back <- data.back.sort[,7]

energies.obs <- energies.obs*(-1)
energies.back <- energies.back*(-1)

M <- fgev(energies.back, std.err = FALSE)

l <- fitted.values(M)[1]
sc <- fitted.values(M)[2]
sh <- fitted.values(M)[3]

intarna.pvalues <- 1-pgev(energies.obs, loc=l, scale=sc, shape=sh)

#intarna.pvalues.corr <- p.adjust(intarna.pvalues, method = "BH")
#intarna.pvalues <- 1-pfrechet(energies.obs.pos, loc=l, scale=sc, shape=sh)

data.obs.sort[,"pvalues"] <- intarna.pvalues

write.table(file=paste(args[1],".pvalues", sep = ""), data.obs.sort, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
