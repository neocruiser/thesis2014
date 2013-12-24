source("http://bioconductor.org/biocLite.R")
biocLite("impute")

?samr

#########################################################################
#								SECTION -1-								#
#########################################################################
## PREPARING OF SUBSET A
## Significance analysis of microarray (SAM)
## Soft clustering afterwards
library(samr)
{
## Imputation for microarray
## SAMPLE DESCRIPTION
y = paste(rep(1,15), "Time", rep(seq(1,5),3), sep="")
start=c(1,6,11)
for(i in start){
y[i]=paste(y[i],"Start",sep="")
}
for(i in  start+4){
y[i]=paste(y[i],"End",sep="")
}
## Preprocessing Genes wheighted w+2
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)	## W/n+Weight
			#setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
			#genes <- read.table("meansW2.txt", sep="\t", header=T)
expressed <- as.matrix(genes[,3:17])
x <- as.matrix(expressed[,c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15)])
## PUT ALL THE data together
## Only 8 permutations are run
## I tried to average the samples. Warning: 1 timecourse in 1 class (unreliable)
set.seed(98231)
data_sam <- list(x=x,y=y, geneid=as.character(genes$X),genenames=as.character(genes$locus), logged2=TRUE)
samr.obj<- samr(data_sam, resp.type="One class timecourse",nperms=20, time.summary.type="slope")
## COmpute the nb of sig genes and the FDRs
delta.table <- samr.compute.delta.table(samr.obj)
cat("\n","\n","Selected SAM FDRs","\n")
round(delta.table[c(24,27,30,32,34,37,39),-c(3,6:8)],3)
delta.table
}
## Choose the value of delta & PLOT
delta=.59
{
samr.plot(samr.obj,delta)
siggenes.table<-samr.compute.siggenes.table(samr.obj, delta, data_sam, delta.table)
## EXTRACT
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
write.table(siggenes.table$genes.up, "samr_siggenes_up.txt", sep="\t", quote=F)
write.table(siggenes.table$genes.lo, "samr_siggenes_lo.txt", sep="\t", quote=F)
## check redundancy
str(siggenes.table)
cat("\n","\n","Positive GENES","\n")
head(siggenes.table$genes.up)
tail(siggenes.table$genes.up)
cat("\n","\n","Negative GENES","\n")
head(siggenes.table$genes.lo)
tail(siggenes.table$genes.lo)
}
## DATASET A
## Extract logs
## Downregulated
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)	## W/n+Weight
quanti_lo <- siggenes.table$genes.lo[,2]
genes_lo <- genes[genes$locus %in% quanti_lo,]
genes_lo <- t(genes_lo[,-c(1:2)])
## Upregulated
quanti_up <- siggenes.table$genes.up[,2]
genes_up <- genes[genes$locus %in% quanti_up,]
genes_up <- t(genes_up[,-c(1:2)])
dim(genes_lo);dim(genes_up)
## Averaging Samples (by col)
## Extract logs
x_lo <- data.frame(Y=paste(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))),genes_lo)
x_up <- data.frame(Y=paste(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))),genes_up)
library(plyr)
y_lo <- ddply(x_lo, .(Y), numcolwise(mean))
y_up <- ddply(x_up, .(Y), numcolwise(mean))
## Final dataset (average col)
average_lo <- t(y_lo[,-1])
average_up <- t(y_up[,-1])
colnames(average_lo) <- c("Egg", "Trocophore", "Veliger", "Pediveliger", "Juvenile")
colnames(average_up) <- c("Egg", "Trocophore", "Veliger", "Pediveliger", "Juvenile")
head(average_lo)
dim(average_lo);dim(average_up)
## Append the name (locus) to rownames
genes_lo <- genes[genes$locus %in% quanti_lo,]
rownames(average_lo) <- genes_lo$locus
genes_up <- genes[genes$locus %in% quanti_up,]
rownames(average_up) <- genes_up$locus
}
## ******************** START ************************
## Soft clustering
## Create ExpressionSet
## Find clusters
{
library(Biobase)
## create a minimal Expression Set
minimalSet_ns_lo <- ExpressionSet(assayData=average_lo)
minimalSet_ns_up <- ExpressionSet(assayData=average_up)
## Soft clustering and visualisation
library(Mfuzz)
## Standardization (relative expression vs absolute ampliture)
minimalSet_lo <- standardise(minimalSet_ns_lo)
minimalSet_up <- standardise(minimalSet_ns_up)
## Fuzzifier m testing
cat("\n\n 'm' estimate for LOW genes",mestimate(minimalSet_lo),"\n")
cat("\n\n 'm' estimate for UP genes",mestimate(minimalSet_up),"\n")
}
## ******************** HOLD *************************
## set c and m before RUN
## RUN soft clustering *************LOW************
cl_lo <- mfuzz(minimalSet_lo,c=3,m=1.5)
mfuzz.plot(minimalSet_lo,cl=cl_lo,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
{
## coupling between clusters (high values between clusters = smilar pattern)
O_lo <- overlap(cl_lo)
over_plot_lo <- overlap.plot(cl_lo, over=O_lo, thres=0.05)
cat("\n","\n","If coupling is large, cluster patterns are similar ---- LOW","\n")
over_plot_lo
## Analysis of global clustering structures
cat("\n","\n","LOW genes","\n")
cl2_lo <- mfuzz(minimalSet_lo,c=3,m=1.5)
mfuzz.plot(minimalSet_lo,cl=cl2_lo,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
O2_lo <- overlap(cl2_lo)
over_plot_lo <- overlap.plot(cl2_lo, over=O2_lo,P=over_plot_lo, thres=0.05)
## Additional
## See available components & extract GENES
cat("\n","\n","LOW genes","\n");cl_lo
## Plot the center of the third cluster
plot(cl_lo[[1]][3,],type="l",ylab="Expression")
## Extract the first 10 membership of the first cluster
#cl_lo[[4]][1:10,1]
## Extract genes forming the alpha cores of soft clusters
cat("\n","\n","LOW genes","\n")
acore_list_lo <- acore(minimalSet_lo, cl=cl_lo, min.acore=0.6)
str(acore_list_lo)
}
## ******************** HOLD *************************
## set c and m before RUN
## RUN soft clustering *************UP************
cl_up <- mfuzz(minimalSet_up,c=15,m=1.5)
mfuzz.plot(minimalSet_up,cl=cl_up,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
{
## coupling between clusters (high values between clusters = smilar pattern)
O_up <- overlap(cl_up)
over_plot_up <- overlap.plot(cl_up, over=O_up, thres=0.05)
cat("\n","\n","If coupling is large, cluster patterns are similar ----- UP","\n")
over_plot_up
## Analysis of global clustering structures
cat("\n","\n","UP genes","\n")
cl2_up <- mfuzz(minimalSet_up,c=15,m=1.5)
mfuzz.plot(minimalSet_up,cl=cl2_up,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
O2_up <- overlap(cl2_up)
over_plot_up <- overlap.plot(cl2_up, over=O2_up,P=over_plot_up, thres=0.05)
## Additional
## See available components & extract GENES
cat("\n","\n","UP genes","\n");cl_up
## Plot the center of the third cluster
plot(cl_up[[1]][3,],type="l",ylab="Expression")
## Extract the first 10 membership of the first cluster
#cl_up[[4]][1:10,1]
## Extract genes forming the alpha cores of soft clusters
cat("\n","\n","UP genes","\n")
acore_list_up <- acore(minimalSet_up, cl=cl_up, min.acore=0.6)
str(acore_list_up)
}
## Extract Locus for BALST2GO search ID
{
setA_lo <- NULL
for(i in 1:length(acore_list_lo)) setA_lo <- rbind(setA_lo, matrix(acore_list_lo[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(setA_lo),"\n")
setA_up <- NULL
for(i in 1:length(acore_list_up)) setA_up <- rbind(setA_up, matrix(acore_list_up[[i]]$NAME))
cat("\n\n number of clustered UP genes",nrow(setA_up),"\n")
}


setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
lsos(pat="set.*")
save(list=ls(pattern="set.*"),file="sets.Rdata")

#########################################################################
#								SECTION -2-								#
#########################################################################
## PATTERN discovery
## SUBSET B
## Imputation for microarray

## Preprocessing Genes wheighted w+2
## ***************HOLD************ Set eigengenes 'y'. (eig1: d=1.32) ; (eig2: d=1.23) ; (eig3: d=1.16)
y <- c(3)
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
expressed <- as.matrix(genes[,3:17])
x <- as.matrix(expressed[,c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15)])
## Set the eigengene. START w/ 1,2 ...
## Compute SAM
library(samr)
set.seed(98234)
data_pat <- list(x=x,eigengene.number=y, geneid=as.character(genes$X),genenames=as.character(genes$locus), logged2=TRUE)
pat_obj<- samr(data_pat, resp.type="Pattern discovery",nperms=100)
delta.table <- samr.compute.delta.table(pat_obj)
cat("\n","\n", "Eigengenes and their FDRs","\n")
delta.table
}
## ********HOLD***************
## SET DELTA
delta=1.12
{
samr.plot(pat_obj,delta)
siggenes.table<-samr.compute.siggenes.table(pat_obj, delta, data_pat, delta.table)
##EXTRACT summary table
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
write.table(siggenes.table$genes.up, "pattern_siggenes_up_eig3.txt", sep="\t", quote=F)
write.table(siggenes.table$genes.lo, "pattern_siggenes_lo_eig3.txt", sep="\t", quote=F)

## check redundancy
str(siggenes.table)
cat("\n","\n","Positive GENES","\n")
head(siggenes.table$genes.up)
tail(siggenes.table$genes.up)
cat("\n","\n","Negative GENES","\n")
head(siggenes.table$genes.lo)
tail(siggenes.table$genes.lo)
}
## DATASET B
## Extract logs for UP and LOW
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
quanti_lo <- siggenes.table$genes.lo[,2]
genes_lo <- genes[genes$locus %in% quanti_lo,]
genes_lo <- t(genes_lo[,-c(1:2)])
## Upregulated
quanti_up <- siggenes.table$genes.up[,2]
genes_up <- genes[genes$locus %in% quanti_up,]
genes_up <- t(genes_up[,-c(1:2)])
dim(genes_lo);dim(genes_up)
## Averaging Samples (by col)
## Extract logs
x_lo <- data.frame(Y=paste(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))),genes_lo)
x_up <- data.frame(Y=paste(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))),genes_up)
library(plyr)
y_lo <- ddply(x_lo, .(Y), numcolwise(mean))
y_up <- ddply(x_up, .(Y), numcolwise(mean))
## Final dataset (average col)
average_lo <- t(y_lo[,-1])
average_up <- t(y_up[,-1])
colnames(average_lo) <- c("Egg", "Trocophore", "Veliger", "Pediveliger", "Juvenile")
colnames(average_up) <- c("Egg", "Trocophore", "Veliger", "Pediveliger", "Juvenile")
head(average_lo)
dim(average_lo);dim(average_up)
## Append the name (locus) to rownames
genes_lo <- genes[genes$locus %in% quanti_lo,]
rownames(average_lo) <- genes_lo$locus
genes_up <- genes[genes$locus %in% quanti_up,]
rownames(average_up) <- genes_up$locus
}


## ******************** START ************************
## Soft clustering
## Create ExpressionSet
## Find clusters
{
library(Biobase)
## create a minimal Expression Set
minimalSet_ns_lo <- ExpressionSet(assayData=average_lo)
minimalSet_ns_up <- ExpressionSet(assayData=average_up)
## Soft clustering and visualisation
library(Mfuzz)
## Standardization (relative expression vs absolute ampliture)
minimalSet_lo <- standardise(minimalSet_ns_lo)
minimalSet_up <- standardise(minimalSet_ns_up)
## Fuzzifier m testing
cat("\n\n 'm' estimate for LOW genes",mestimate(minimalSet_lo),"\n")
cat("\n\n 'm' estimate for UP genes",mestimate(minimalSet_up),"\n")
}
## ******************** HOLD *************************
## set c and m before RUN
## RUN soft clustering *************LOW************
cl_lo <- mfuzz(minimalSet_lo,c=5,m=1.5)
mfuzz.plot(minimalSet_lo,cl=cl_lo,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
{
## coupling between clusters (high values between clusters = smilar pattern)
O_lo <- overlap(cl_lo)
over_plot_lo <- overlap.plot(cl_lo, over=O_lo, thres=0.05)
cat("\n","\n","If coupling is large, cluster patterns are similar ---- LOW","\n")
over_plot_lo
## Analysis of global clustering structures
cat("\n","\n","LOW genes","\n")
cl2_lo <- mfuzz(minimalSet_lo,c=5,m=1.5)
mfuzz.plot(minimalSet_lo,cl=cl2_lo,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
O2_lo <- overlap(cl2_lo)
over_plot_lo <- overlap.plot(cl2_lo, over=O2_lo,P=over_plot_lo, thres=0.05)
## Additional
## See available components & extract GENES
cat("\n","\n","LOW genes","\n");cl_lo
## Plot the center of the third cluster
plot(cl_lo[[1]][3,],type="l",ylab="Expression")
## Extract the first 10 membership of the first cluster
#cl_lo[[4]][1:10,1]
## Extract genes forming the alpha cores of soft clusters
cat("\n","\n","LOW genes","\n")
acore_list_lo <- acore(minimalSet_lo, cl=cl_lo, min.acore=0.6)
str(acore_list_lo)
}
## ******************** HOLD *************************
## set c and m before RUN
## RUN soft clustering *************UP************
cl_up <- mfuzz(minimalSet_up,c=9,m=1.5)
mfuzz.plot(minimalSet_up,cl=cl_up,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
{
## coupling between clusters (high values between clusters = smilar pattern)
O_up <- overlap(cl_up)
over_plot_up <- overlap.plot(cl_up, over=O_up, thres=0.05)
cat("\n","\n","If coupling is large, cluster patterns are similar ----- UP","\n")
over_plot_up
## Analysis of global clustering structures
cat("\n","\n","UP genes","\n")
cl2_up <- mfuzz(minimalSet_up,c=9,m=1.5)
mfuzz.plot(minimalSet_up,cl=cl2_up,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
O2_up <- overlap(cl2_up)
over_plot_up <- overlap.plot(cl2_up, over=O2_up,P=over_plot_up, thres=0.05)
## Additional
## See available components & extract GENES
cat("\n","\n","UP genes","\n");cl_up
## Plot the center of the third cluster
plot(cl_up[[1]][3,],type="l",ylab="Expression")
## Extract the first 10 membership of the first cluster
cl_up[[4]][1:10,1]
## Extract genes forming the alpha cores of soft clusters
cat("\n","\n","UP genes","\n")
acore_list_up <- acore(minimalSet_up, cl=cl_up, min.acore=0.7)
str(acore_list_up)
}
## Extract Locus for BALST2GO search ID
## eigengene 1
{
setB_lo_eig1 <- NULL
for(i in 1:length(acore_list_lo)) setB_lo_eig1 <- rbind(setB_lo_eig1, matrix(acore_list_lo[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(setB_lo_eig1),"\n")
setB_up_eig1 <- NULL
for(i in 1:length(acore_list_up)) setB_up_eig1 <- rbind(setB_up_eig1, matrix(acore_list_up[[i]]$NAME))
cat("\n\n number of clustered UP genes",nrow(setB_up_eig1),"\n")
}
## eigengene 2
{
setB_lo_eig2 <- NULL
for(i in 1:length(acore_list_lo)) setB_lo_eig2 <- rbind(setB_lo_eig2, matrix(acore_list_lo[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(setB_lo_eig2),"\n")
setB_up_eig2 <- NULL
for(i in 1:length(acore_list_up)) setB_up_eig2 <- rbind(setB_up_eig2, matrix(acore_list_up[[i]]$NAME))
cat("\n\n number of clustered UP genes",nrow(setB_up_eig2),"\n")
}
## eigengene 3
{
setB_lo_eig3 <- NULL
for(i in 1:length(acore_list_lo)) setB_lo_eig3 <- rbind(setB_lo_eig3, matrix(acore_list_lo[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(setB_lo_eig3),"\n")
setB_up_eig3 <- NULL
for(i in 1:length(acore_list_up)) setB_up_eig3 <- rbind(setB_up_eig3, matrix(acore_list_up[[i]]$NAME))
cat("\n\n number of clustered UP genes",nrow(setB_up_eig3),"\n")
}

setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
lsos(pat="set.*")
save(list=ls(pattern="set.*"),file="sets.Rdata")

## OPTIONAL
##plot x=called; y=FDR median
## CHANGE the x_lim and y_lim on every eigengene
{
x_lim <- subset(delta.table, delta.table[,4]<2100, select="# called")
y_lim <- subset(delta.table, delta.table[,2]<70, select="# med false pos")
plot(x_lim, y_lim, type="s", las=1, xlab="log ratio #called genes", ylab="FDR", frame=F, log="x")
axis(side=1, at=c(10,2000), labels=F)
axis(side=2, at=c(0, 100),labels=F)
}

#########################################################################
#								SECTION -3-								#
#########################################################################
## Specific PATTERN DISCOVERY using the QUANTITATIVE option
## SUBSET C
## SOft clustering Afterwards
## Imputation for microarray

## SAMPLE DESCRIPTION for : (E T) (V P) J
## Already tested : (21000 d=0.4 133);(11000 d=0.31 27 5);(32100 d=0.49 318) = SETC1. (00012 d=0.68);(00001 d=063)=SETC2
y <- paste(rep(c(0,0,0,0,1),3), sep="")
{
library(samr)
## Preprocessing Genes wheighted w+2
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
expressed <- as.matrix(genes[,3:17])
x <- as.matrix(expressed[,c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15)])
## Compute SAM
set.seed(1325)
data_pat <- list(x=x,y=y, geneid=as.character(genes$X),genenames=as.character(genes$locus), logged2=TRUE)
pat_obj<- samr(data_pat, resp.type="Quantitative",nperms=50, regression.method="ranks")
delta.table <- samr.compute.delta.table(pat_obj)
cat("\n","\n", "Eigengenes and their FDRs","\n")
delta.table
}
##**********HOLD*************
## SET DELTA
delta=0.58
{
#samr.plot(pat_obj,delta)
siggenes.table<-samr.compute.siggenes.table(pat_obj, delta, data_pat, delta.table)
## check redundancy
str(siggenes.table)
cat("\n\nPositive GENES\n")
head(siggenes.table$genes.up)
tail(siggenes.table$genes.up)
cat("\n","\n","Negative GENES","\n")
head(siggenes.table$genes.lo)
tail(siggenes.table$genes.lo)
}
## Add suffix to create a file for every pattern
quanti_siggenes_up21000 <- siggenes.table$genes.up
quanti_siggenes_lo21000 <- siggenes.table$genes.lo
quanti_siggenes_up11000 <- siggenes.table$genes.up
quanti_siggenes_lo11000 <- siggenes.table$genes.lo
quanti_siggenes_up32100 <- siggenes.table$genes.up
quanti_siggenes_lo32100 <- siggenes.table$genes.lo
quanti_siggenes_up00012 <- siggenes.table$genes.up
quanti_siggenes_lo00012 <- siggenes.table$genes.lo
quanti_siggenes_up00001 <- siggenes.table$genes.up
quanti_siggenes_lo00001 <- siggenes.table$genes.lo

setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
lsos(pat="quanti_sigg.*|set.*")
save(list=ls(pattern="quanti_sigg.*|set.*"),file="sets.Rdata")

## *************STOP***************
## PREPARE CLUSTERING 
## Extract logs **********LOW****************
### set variaX to one of the previous saved sets 
variaX <- quanti_siggenes_up00001[,2]
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
genes_lo <- genes[genes$locus %in% variaX,]
genes_lo <- t(genes_lo[,-c(1:2)])
#x <- colnames(genes_lo[,4:18])
## Plot all lines
plot(genes_lo[,1], type="l", xlab="Developmental stages", ylab="log ratio")
for(j in 2:ncol(genes_lo)) lines(genes_lo[,j], type="l", col=j, lwd=0.5)
## Averaging Samples (by col)
## Extract logs
x <- data.frame(Y=paste(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))),genes_lo)
library(plyr)
y <- ddply(x, .(Y), numcolwise(mean))
## Final dataset (average col)
average_lo <- t(y[,-1])
colnames(average_lo) <- c("Egg", "Trocophore", "Veliger", "Pediveliger", "Juvenile")
genes_lo <- genes[genes$locus %in% variaX,]
rownames(average_lo) <- genes_lo$locus
head(average_lo)
}
## create a minimal ExpressionSet & run MFUZZ
{
library(Biobase)
minimalSet_ns <- ExpressionSet(assayData=average_lo)
## Soft clustering and visualisation
library(Mfuzz)
## Standardization (relative expression vs absolute ampliture)
minimalSet <- standardise(minimalSet_ns)
## Fuzzifier m testing
mestimate(minimalSet)
## RUN soft clustering
cl <- mfuzz(minimalSet,c=12,m=1.5)
mfuzz.plot(minimalSet,cl=cl,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
}
## coupling between clusters (high values between clusters = smilar pattern)
{
O <- overlap(cl)
over_plot <- overlap.plot(cl, over=O, thres=0.05)
cat("\n","\n","If coupling is large, cluster patterns are similar","\n")
over_plot
## Analysis of global clustering structures
cl2 <- mfuzz(minimalSet,c=12,m=1.5)
mfuzz.plot(minimalSet,cl=cl2,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
O2 <- overlap(cl2)
over_plot <- overlap.plot(cl2, over=O2,P=over_plot, thres=0.05)
}
## ADDITIONAL OPTIONS (EXTRACT CLUSTERS)
## See available components & extract GENES
{
cl
## Plot the center of the third cluser
plot(cl[[1]][3,],type="l",ylab="Expression")
## Extract the first 10 membership of the first cluster
cl[[4]][1:10,1]
## Extract genes forming the alpha cores of soft clusters
acore_list <- acore(minimalSet, cl=cl, min.acore=0.6)
str(acore_list)
}

## (21000)
{
setC_up_21000 <- NULL
for(i in 1:length(acore_list)) setC_up_21000 <- rbind(setC_up_21000, matrix(acore_list[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(setC_up_21000),"\n")
}
## (21000)
{
setC_lo_21000 <- NULL
for(i in 1:length(acore_list)) setC_lo_21000 <- rbind(setC_lo_21000, matrix(acore_list[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(setC_lo_21000),"\n")
}
## (11000)
{
setC_up_11000 <- NULL
for(i in 1:length(acore_list)) setC_up_11000 <- rbind(setC_up_11000, matrix(acore_list[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(setC_up_11000),"\n")
}
## (11000)
{
setC_lo_11000 <- NULL
for(i in 1:length(acore_list)) setC_lo_11000 <- rbind(setC_lo_11000, matrix(acore_list[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(setC_lo_11000),"\n")
}
## (32100)
{
setC_lo_32100 <- NULL
for(i in 1:length(acore_list)) setC_lo_32100 <- rbind(setC_lo_32100, matrix(acore_list[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(setC_lo_32100),"\n")
}
## (00012)
{
setC_up_00012 <- NULL
for(i in 1:length(acore_list)) setC_up_00012 <- rbind(setC_up_00012, matrix(acore_list[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(setC_up_00012),"\n")
}
## (00001)
{
setC_up_00001 <- NULL
for(i in 1:length(acore_list)) setC_up_00001 <- rbind(setC_up_00001, matrix(acore_list[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(setC_up_00001),"\n")
}

setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
lsos(pat="quanti_sigg.*|set.*")
save(list=ls(pattern="quanti_sigg.*|set.*"),file="sets.Rdata")


## Extract unique Subsets C, for paper 2
## CombineTemplate patterns
setC <- rbind(setC_lo_11000,setC_lo_21000,setC_lo_32100,setC_up_00001,setC_up_00012,setC_up_21000)
subsetC <- as.matrix(setC[!duplicated(setC[,1]),])
## clustering of subset C
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
genes_lo <- genes[genes$locus %in% subsetC[,1],]
genes_lo <- t(genes_lo[,-c(1:2)])
#x <- colnames(genes_lo[,4:18])
## Plot all lines
plot(genes_lo[,1], type="l", xlab="Developmental stages", ylab="log ratio")
for(j in 1:ncol(genes_lo)) lines(genes_lo[,j], type="l", col=j, lwd=0.5)
## Averaging Samples (by col)
## Extract logs
x <- data.frame(Y=paste(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))),genes_lo)
library(plyr)
y <- ddply(x, .(Y), numcolwise(mean))
## Final dataset (average col)
average <- t(y[,-1])
colnames(average) <- c("Egg", "Trocophore", "Veliger", "Pediveliger", "Juvenile")
genes_lo <- genes[genes$locus %in% subsetC[,1],]
rownames(average) <- genes_lo$locus
head(average)
}
## create a minimal ExpressionSet & run MFUZZ
## ********************* subset C *****************
{
library(Biobase)
minimalSet_ns <- ExpressionSet(assayData=average)
## Soft clustering and visualisation
library(Mfuzz)
## Standardization (relative expression vs absolute ampliture)
minimalSet <- standardise(minimalSet_ns)
## Fuzzifier m testing
mestimate(minimalSet)
}
## RUN soft clustering
cl <- mfuzz(minimalSet,c=20,m=1.5)
mfuzz.plot(minimalSet,cl=cl,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
## coupling between clusters (high values between clusters = smilar pattern)
## ************************** Subset C**********************
{
O <- overlap(cl)
over_plot <- overlap.plot(cl, over=O, thres=0.05)
cat("\n","\n","If coupling is large, cluster patterns are similar","\n")
over_plot
## Analysis of global clustering structures
cl2 <- mfuzz(minimalSet,c=20,m=1.5)
mfuzz.plot(minimalSet,cl=cl2,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
O2 <- overlap(cl2)
over_plot <- overlap.plot(cl2, over=O2,P=over_plot, thres=0.05)
}
## ADDITIONAL OPTIONS (EXTRACT CLUSTERS)
## See available components & extract GENES
## **************Subset C
{
cl
## Plot the center of the third cluser
plot(cl[[1]][3,],type="l",ylab="Expression")
## Extract the first 10 membership of the first cluster
cl[[4]][1:10,1]
## Extract genes forming the alpha cores of soft clusters
acore_list <- acore(minimalSet, cl=cl, min.acore=0.6)
str(acore_list)


}
## Extract Locus for BALST2GO search ID
## SubsetC final
{
subsetC_final <- NULL
for(i in 1:length(acore_list)) subsetC_final <- rbind(subsetC_final, matrix(acore_list[[i]]$NAME))
cat("\n\n number of clustered LOW genes",nrow(subsetC_final),"\n")
}

setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
lsos(pat="quanti_sigg.*|set.*")
save(list=ls(pattern="quanti_sigg.*|set.*"),file="sets.Rdata")

#########################################################################
#								SECTION -4-								#
#########################################################################
## clustering of subsets of all previous sections

{					## Load SoftClustering Mfuzz package --- 12 clusters
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
load("clusters12_Lm.Rdata", .GlobalEnv)
library(Mfuzz)
lsos(pat="")
}
set.seed(349815720)

# LOAD ALL SUBSETS A B C
{
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
load("sets.Rdata", .GlobalEnv)
lsos(pat="quanti_sigg.*|set.*")
}

## merge all sets into one file, then load it BALST2GO
## LIBRARY 1
all.set <- as.matrix(rbind(setA_lo,setA_up,setB_lo_eig1,setB_lo_eig2,setB_lo_eig3,setB_up_eig1,setB_up_eig2,setB_up_eig3,subsetC_final))
nrow(all.set)
all.set <- all.set[!duplicated(all.set[,1]),]
length(all.set)
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
write.table(all.set,"allSet_Lm.txt",sep="\t",quote=F)
## Remove additional entries in NOTEPAD++ ********************IMPORTANT****************
## RUN BLAST2GO
## LOAD Medulsi project
## load all.sets
## Launch a Fishers enrishment analysis
## copy paste 1 GO set (Test group) into .txt, then replace the: , with: \r\n. Then read in the file

## OPTIONAL Load all GO sets (if you want to cluster them)
{
GO1 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0043232.txt", sep="\t", header=F)
GO2 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0006412.txt", sep="\t", header=F)
GO3 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0010467.txt", sep="\t", header=F)
GO4 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0030529.txt", sep="\t", header=F)
GO5 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0005694.txt", sep="\t", header=F)
GO6 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0005198.txt", sep="\t", header=F)
GO7 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0032991.txt", sep="\t", header=F)
GO8 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0044424.txt", sep="\t", header=F)
lsos(pat="GO.")
}
## Extract logs
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
all.set <- read.table("allSet_Lm.txt", sep="\t", header=T)
## SET k
k=all.set
{					# Extract data
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
genes.c <- genes[genes$locus %in% k[,1],]
genes.c <- t(genes.c[,-c(1:2)])
## Averaging Samples (by col)
## Extract logs
x <- data.frame(Y=paste(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))),genes.c)
library(plyr)
y <- ddply(x, .(Y), numcolwise(mean))
## Final dataset (average col)
average <- t(y[,-1])
colnames(average) <- c("Egg", "Trocophore", "Veliger", "Pediveliger", "Juvenile")
genes.c <- genes[genes$locus %in% k[,1],]
rownames(average) <- genes.c$locus
head(average)
}
## create a minimal ExpressionSet & run MFUZZ
## ********************* GO Term *****************
## SET i and j
i=12 ; j=1.55
{
library(Biobase)
minimalSet_ns <- ExpressionSet(assayData=average)
## Soft clustering and visualisation
library(Mfuzz)
## Standardization (relative expression vs absolute ampliture)
minimalSet <- standardise(minimalSet_ns)
## Fuzzifier m testing
cat("\n\n 'm' estimate for genes",mestimate(minimalSet),"\n")
## RUN soft clustering
cl <- mfuzz(minimalSet,c=i,m=j)
mfuzz.plot(minimalSet,cl=cl,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
}
## coupling between clusters (high values between clusters = smilar pattern)
{
O <- overlap(cl)
over_plot <- overlap.plot(cl, over=O, thres=0.05)
cat("\n","\n","If coupling is large, cluster patterns are similar","\n")
over_plot
## Analysis of global clustering structures
cl2 <- mfuzz(minimalSet,c=i,m=j)
mfuzz.plot(minimalSet,cl=cl2,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
O2 <- overlap(cl2)
over_plot <- overlap.plot(cl2, over=O2,P=over_plot, thres=0.05)
}
## ADDITIONAL OPTIONS (EXTRACT CLUSTERS)
## See available components & extract GENES

## Extract genes forming the alpha cores of soft clusters
acore_list09<- acore(minimalSet, cl=cl, min.acore=0.9)
str(acore_list07)
for(i in 1:length(acore_list07)) print(nrow(acore_list07[[i]]))

## Extract membership
{
cluster1 <- data.frame(x=acore_list07[[1]]$NAME)
cluster2 <- data.frame(x=acore_list07[[2]]$NAME)
cluster3 <- data.frame(x=acore_list07[[3]]$NAME)
cluster4 <- data.frame(x=acore_list07[[4]]$NAME)
cluster5 <- data.frame(x=acore_list07[[5]]$NAME)
cluster6 <- data.frame(x=acore_list07[[6]]$NAME)
cluster7 <- data.frame(x=acore_list07[[7]]$NAME)
cluster8 <- data.frame(x=acore_list07[[8]]$NAME)
cluster9 <- data.frame(x=acore_list07[[9]]$NAME)
cluster10 <- data.frame(x=acore_list07[[10]]$NAME)
cluster11 <- data.frame(x=acore_list07[[11]]$NAME)
cluster12 <- data.frame(x=acore_list07[[12]]$NAME)
}

## Plot the center of the clusters
par(mfrow=c(3,3))
********* History > Recording (after flashing par())
for(z in 1:i) plot(cl[[1]][z,],type="l",ylab="Expression")


## Save a copy of the memory
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
lsos(pat="cluster.*|acore_list.*")
save(list=ls(pattern="cluster.*|acore_list.*"),file="clusters12_Lm.Rdata")




