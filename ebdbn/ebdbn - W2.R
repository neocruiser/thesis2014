## Rau 2009
## EBDBN

#########################################################################
#								SECTION -1-								#
#########################################################################
All.set subset A, B, C (SAM)
==================================Load(0)

{					## Load L.Rdata
setwd("C:/Dropbox/Workshop2013/Work/R/ebdbn/")
load("networks.Rdata", .GlobalEnv)
require(ebdbNet)
library(GeneNet)
require(preprocessCore)
tmp <- runif(1) ## Initialize random number generator
set.seed(45684598) ## Set seed
lsos(pat="")
}
==================================RUN(a)(b)(c)
## Transcription factors (Molecular function)
## Load Transcription factors TF.Rdata extracted from Library 1 (after molecular GO combined graph in Blast2Go)
{
setwd("C:/Dropbox/Workshop2013/Work/R/ebdbn/")
load("TF.Rdata", .GlobalEnv)
print(lsos(pat="GO.*"))

helloTF <- rbind(GO_0003700,GO_0030234,GO_0030528)
dim(helloTF)
helloTF <- helloTF[! duplicated(helloTF),]
GOx3 <- data.frame(V1=helloTF)
cat("\nDimentions of GOx3 rbind 3 GO sets before",dim(GOx3),"\n")
}
## Choose GO
TFGO=GOx3
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
genes_SIG <- genes_coc[genes_coc$locus %in% TFGO$V1,]
genes_mat <- as.matrix(genes_SIG[,4:18])

## Transform & Normalize
Tgo <- normalize.quantiles(genes_mat)		# IMP: matrix where columns are arrays, rows are probes

rownames(Tgo) <- paste(genes_SIG$locus,"|",genes_SIG$description)
Tgo <- t(Tgo)
cat("Dimentions of the Transcription Factor Matrix :",dim(Tgo),"\n")	#15x1444

speciesF <- as.longitudinal(Tgo, repeats=3, time=c(0,6,12,22,32))
TF <- dataFormat(speciesF)		# Transcription factors
}
==================================
## ALL.SET
## all.set for subset A, B, C (OPTIONAL)
{
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
all.set <- read.table("all.set.txt", sep="\t", header=F)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
genes_SIG <- genes_coc[genes_coc$locus %in% all.set$V2,]
## Remove TF (transcription factors) from all.set
Nsig <- genes_SIG[! genes_SIG[,3] %in% TFGO$V1,]
N <- as.matrix(Nsig[,4:18])
cat("Dimentions of the Observed Genes Matrix :",dim(N),"\n")	#15x1444

## Transform & Normalize
Ta <- normalize.quantiles(N) 		# IMP: matrix where columns are arrays, rows are probes

rownames(Ta) <- Nsig$locus
Ta <- t(Ta)

species <- as.longitudinal(Ta, repeats=3, time=c(0,6,12,22,32))
G <- dataFormat(species)		# Observed Genes

}
==================================RUN(a)
## CLUSTERS (12 clusters)
## Load all clusters 12 from (last soft clustering for paper 2 w all.set -- Library 1) (optional)
{					## Load SoftClustering Mfuzz package --- 12 clusters
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
load("clusters12.Rdata", .GlobalEnv)
lsos(pat="cluster.*")
}
## Clusters (should load unique GO first)
N=cluster1
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
N <- genes_coc[genes_coc$locus %in% N$x,]
cat("Before cleaning :", dim(N),"\n")
## Remove TF (transcription factors) from all.set
Nsig <- N[! N[,3] %in% TFGO$V1,]
cat("After cleaning :", dim(Nsig),"\n")

N <- as.matrix(Nsig[,4:18])
cat("Dimentions of the Observed Genes Matrix :",dim(N),"\n")	#15x1444

## Normalize
## Transform & Normalize
Tcl <- normalize.quantiles(N)		# IMP: matrix where columns are arrays, rows are probes

rownames(Tcl) <- Nsig$locus
Tcl <- t(Tcl)
## Transform to class(longitudinal)
## transform for observed genes
species <- as.longitudinal(Tcl, repeats=3, time=c(0,6,12,22,32))
G <- dataFormat(species)		# Observed Genes
}
==================================RUN(a)
## SEARCH
## Search and locate where TF (extracted from Library 1) in GOx3 are in clusters of Library 1 (12 clusters)
{						## RUN (1)
clusters <- list(cluster1,cluster2,cluster3,cluster4,cluster5,cluster6,cluster7,cluster8,cluster9,cluster10,cluster11,cluster12)
Lo=as.vector(as.matrix(GOx3))
Lo
}
## SET (LoX)
LoX=c("Locus554")
{						## RUN (2)
hello <- list()
for(i in 1:length(clusters)) {
hello[[i]] <- grepl(LoX,clusters[[i]]$x) ; 
hello[[i]] <- length(which(hello[[i]])) 
}
str(hello)
}
## (optional) only if TF not found ^ (lower alpha core =0.3)
{						## RUN (2) \\ change the suffix to _03, _05 or _06
hello <- list()
for(i in 1:length(acore_list01)) {
hello[[i]] <- grepl(LoX,acore_list01[[i]]$NAME) ; 
hello[[i]] <- length(which(hello[[i]])) 
}
str(hello)
}
==================================RUN(b)
## Anova one way results
## LOAD 5 subsets (from 5 different Templates)
{
setwd("C:/Dropbox/Workshop2013/Work/R/ANOVA/")
anova_template_EGG <- read.table("template-logs_EGG.txt", sep="\t", header=T)
anova_template_TROCO <- read.table("template-logs_TROCO.txt", sep="\t", header=T)
anova_template_VELI <- read.table("template-logs_VELI.txt", sep="\t", header=T)
anova_template_PEDI <- read.table("template-logs_PEDI.txt", sep="\t", header=T)
anova_template_JUV <- read.table("template-logs_JUV.txt", sep="\t", header=T)
lsos(pat="anova_template_.*")
}
## choose subset
N=anova_template_EGG
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansW2.txt", sep="\t", header=T, row.names=3)
description <- genes[rownames(genes) %in% colnames(N),]
description <- description[ ! rownames(description) %in% TFGO$V1,]
cat("Dimentions of the Description dataframe :",dim(description),"\n")	# mx15

NNmatrix <- as.matrix(t(N))
cat("Dimentions of the Observed Genes Before GO :",dim(NNmatrix),"\n")	# mx15

NNsig <- NNmatrix[! rownames(NNmatrix) %in% TFGO$V1,]
cat("Dimentions of the Observed Genes After GO :",dim(NNsig),"\n")	# mx15

## Transform & Normalize
NNnormalize <- normalize.quantiles(NNsig) 		# IMP: matrix where columns are arrays, rows are probes
	# rownames(NNnormalize) <- rownames(NNsig)
rownames(NNnormalize) <- paste(rownames(description),"|",description$description)
NNnormalize <- t(NNnormalize)
species <- as.longitudinal(NNnormalize, repeats=3, time=c(0,6,12,22,32))
G <- dataFormat(species)		# Observed Genes
}
==================================RUN(c)
## SubsetE (Growth, Development, Fatty acids metabolism --- chosen from Library1, after Process GO combined graph with Blast2Go)
## Load setE
{
setwd("C:/Dropbox/Workshop2013/Work/R/ebdbn/")
load("setE.Rdata", .GlobalEnv)
lsos(pat="setE*")

fooE <- rbind(setEgrowth,setElipid,setEdevelopment)
dim(fooE)
fooE <- fooE[! duplicated(fooE),]
setEx3 <- data.frame(V1=fooE)
cat("\nDimentions of setEx3 rbind 3 Processes sets before",dim(setEx3),"\n")
}
## Choose GO
x3setE=setEx3
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansW2.txt", sep="\t", header=T, row.names=3)
description1 <- genes[rownames(genes) %in% x3setE$V1,]
cat("Dimentions of the Observed Genes Before TF extraction :",dim(description1),"\n")	# mx15
description <- description1[ ! rownames(description1) %in% TFGO$V1,]
cat("Dimentions of the Observed Genes After TF extraction :",dim(description),"\n")	# mx15
description2 <- as.matrix(description[,-c(1,2)])
## Transform & Normalize
NNnormalize <- normalize.quantiles(description2) 		# IMP: matrix where columns are arrays, rows are probes
	# rownames(NNnormalize) <- rownames(NNsig)
rownames(NNnormalize) <- paste(rownames(description),"|",description[,2])
NNnormalize <- t(NNnormalize)
species <- as.longitudinal(NNnormalize, repeats=3, time=c(0,6,12,22,32))
G <- dataFormat(species)		# Observed Genes
}
==================================
## START
## Calculate Hankel Hidden States at Lag=1, 2, 3
K <- hankel(G,lag=1,cutoff=0.9,type="median")$dim ; K

## Feedback networks with/without hidden states
net30 <- ebdbn(input="Feedback",G,K=3,conv.1=0.1, conv.2=0.05, conv.3=0.001, verbose=T)
## Input networks with/without hidden state
net30 <- ebdbn(input=TF,G,K=3, conv.1=0.1, conv.2=0.05, conv.3=0.001, verbose=T)
plot(net30, sig.level=0.99, interactive=T, clarify=T) ; gc() ; tmp <- runif(1) ; set.seed(45684598)
## END

## Extract Positif or negatif regulations (for feedback and input alike)
names(net14$z[1,])	## list gene names (leave it 1) then set the net** number
f=26	## choose number from the list
{			## choose f first then change net** inside the loop
as.matrix(net14$z[,f])		## All genes : fth regulator
#	names(net14$z[,1])	## list gene names (leave it 1)
#	net14$z[1,]		## All regulators : First gene
}

## LOOK if target gene is enhanced or repressed at the stage of network inference
f="Locus554"
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
MM <- read.table("meansW2.txt", sep="\t", header=T)
samples15 <- MM[MM$locus %in% c(f),]
require(plyr)
samples <- as.numeric(gl(5,3,15, label=seq(1:5)))
NNdf  <- data.frame(samples=samples, t(samples15[,-c(1:3)]))
NNm <- ddply(NNdf, .(samples), numcolwise(mean))
NNts <- t(NNm[,-1])
colnames(NNts) <-  c("Egg", "Troco", "Veli", "Pedi", "Juv")
row.names(NNts) <- paste(samples15[,3],"-", samples15[,2])
NNts
}

## SAVED NETWORKS AND DESCRIPTIONS
## for allset dataset: (1444 seq)
net : input | vs GO18 | 01 005 0001 | k=4
net2 : input | vs GO18 | 005 001 001 | k=4
## for anova all (107 seq)
net5 : feedback | EGG | minus GO18 |015 001 001 | K=3
net8 : feedback | PEDI | minus GO18 | 0001 | k=3 | .99

net13 : feedback | JUV | vs GO18 | 01 001 0001 | k=3 | .99

*net6 : feedback | EGG | vs GOx3 | 0001 | k=3 | .99
*net7 : feedback | VELI | vs GOx3 | 01 005 0001 | k=3 | .95
*net21 : feedback | setE | vs GOx3 | 0001 | k=1 | .99

**net14 : input | EGG | vs GOx3 | 015 005 001 | k=3 | .99
**net15 : input | VELI | vs GOx3 | 015 001 001 | k=3 | .99
**net16 : input | JUV | vs GOx3 | 01 005 0001 | k=3 | .99
**net20 : input | setE | vs GOx3 | 01 005 0001 | k=3 | .99
*net30 : feedback | setE | 01 005 0001 | k=3 | .98


## Save a copy of the memory
setwd("C:/Dropbox/Workshop2013/Work/R/ebdbn/")
lsos(pat="net.*")
save(list=ls(pattern="net.*"),file="networks.Rdata")

## Create GO molecular function data frames
setwd("C:/Dropbox/Workshop2013/Work/R/ebdbn/")
setDdevelopment <- read.table("test.txt", sep="\t", header=F)
lsos(pat="setE.*")
setwd("C:/Dropbox/Workshop2013/Work/R/ebdbn/")
save(list=ls(pattern="setE.*"),file="setE.Rdata")

#########################################################################
#								SECTION -2-								#
#########################################################################
Example from manual
ebdbNet 1.2.2 (run on R 3.0.1)
====================

require(ebdbNet)
library(GeneNet)
tmp <- runif(1) ## Initialize random number generator
set.seed(4568818) ## Set seed
## Load T-cell activation data
data(tcell)
tc44 <- combine.longitudinal(tcell.10, tcell.34)
get.time.repeats(tcell.34)
## Put data into correct format for algorithm
## (List, with one matrix per replicate (P rows and T columns)
tcell.dat <- dataFormat(tc44)
## Use only subset of T-cell data for faster example
R <- 20 ## 20 replicates
P <- 10 ## 10 genes
{
tcell.sub.dat <- vector("list", R)
rep.sample <- sample(1:44, R)
for(r in 1:R) {
tcell.sub.dat[[r]] <- tcell.dat[[rep.sample[r]]][sample(1:58, P),]
}
## Choose alternative value of K using hankel if hidden states to be estimated
## K <- hankel(tcell.sub.dat, lag = 1)$dim
## Run algorithm (feedback network, no hidden states)
net <- ebdbn(y = tcell.sub.dat, K = 0, input = "feedback", conv.1 = 0.01,
conv.2 = 0.01, conv.3 = 0.001, verbose = TRUE)
}
## Visualize results: in this example, mostly feedback loops
plot(net, sig.level = 0.5)

#########################################################################
#								SECTION -3-								#
#########################################################################

Hello Sleiman,
I took the Cytoscape formatting capabilities out of the package once I started using igraph to plot the graphs internally in R. Below you can find the function from the earlier version that formats the results of ebdbNet for input into Cytoscape. You might have to play around a bit to fit in with the results of the latest version of ebdbNet, notably zscores should be the $z value from an ebdbn object, and sig should be a matrix of 0's and 1's (of the same dimension as zscore) indicating whether zscore is significant at whatever threshold you want. Let me know if it's not clear or if you have other questions.
Best, Andrea

visualize <- function (zscores, sig, type = "feedback")  {
    row <- dim(zscores)[1]
    col <- dim(zscores)[2]
    network <- matrix(0, nrow = row, ncol = col)
    network <- sign(zscores) * sig
    pos <- matrix(which(network == 1, arr.ind = TRUE)[, 2:1], ncol = 2)
    neg <- matrix(which(network == -1, arr.ind = TRUE)[, 2:1], ncol = 2)
    pos.edges <- cbind(pos, rep(1, length(which(network == 1))))
    neg.edges <- cbind(neg, rep(-1, length(which(network ==  -1))))
    edges <- rbind(pos.edges, neg.edges)
    if(type == "input") {
        edges <- cbind(paste("TF", edges[,1], sep = ""), paste("G", edges[,2], sep = ""), edges[,3])
    }
    return(edges)
}

