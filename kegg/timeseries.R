## Extract Logs ; 3 replicates per stage
## set the name and path of the locus containing file
{
setwd("C:/Dropbox/Workshop2013/Work/R/kegg/")
x= read.table("aatRNA.txt", sep="\t", header=F)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
NN <- genes_coc[genes_coc$locus %in% x$V1,]
cat("Dimentions of the Observed Genes Matrix :",dim(NN),"\n")

NNp <- as.matrix(NN[,4:18]) ; rownames(NNp) <- NN$locus
cat("Dimentions of the Observed Genes Matrix :",dim(NNp),"\n")

## Logs with mean per stage
library(plyr)
samples <- as.numeric(gl(5,3,15, label=seq(1:5)))
NNdf  <- data.frame(samples=samples, t(NNp))
NNm <- ddply(NNdf, .(samples), numcolwise(mean))
NNts <- t(NNm[,-1])
cat("Dimentions of the Observed Genes Matrix :",dim(NNts),"\n")
}

y <- rainbow(32)
plot(NNts[1,], type="n",lwd=2, ylim=c(-1.5,1.5))
for(i in 2:nrow(NNts)) lines(NNts[i,], type="l", lwd=2, col=y[i])