


#########################################################################
#								SECTION -1-								#
#########################################################################

library(Biobase)
vignette(package="Biobase")
vignette(package="Biobase", "ExpressionSetIntroduction")


## Preprocessing Genes wheighted w+2
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansW2.txt", sep="\t", header=T)
expressed <- as.matrix(genes[,4:18])
rownames(expressed) <- genes$locus
## arrange the sample columns
#exprs <- as.matrix(expressed[,c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15)])

## Average the sample columns
x <- data.frame(Y=paste(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))),t(expressed))
library(plyr)
y <- ddply(x, .(Y), numcolwise(mean))

## Final dataset (average col)
exprs_m <- t(y[,-1])



## create a minimal Expression Set
minimalSet <- ExpressionSet(assayData=exprs_m)
