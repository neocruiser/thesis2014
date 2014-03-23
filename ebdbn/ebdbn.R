## Rau 2009
## EBDBN

#########################################################################
#								SECTION -1-								#
#########################################################################
## BAYESIAN NETWORKS INFERRED FRO SUBSET D, SUBSET E, LIBRARY 1 WITH REGULATORS (TFS)
==================================Load(always)

{					## Load L.Rdata
setwd("C:/Dropbox/Workshop2013/Work/R/ebdbn/")
load("Networks2.3.Rdata", .GlobalEnv)
require(ebdbNet)
library(GeneNet)
require(preprocessCore)
lsos(pat="network.*")
}
tmp <- runif(1); set.seed(426924)
==================================RUN(a)(b)(c)
## Transcription factors (Molecular function)
## Load Transcription factors TF.Rdata extracted from Library 1 (after molecular GO combined graph in Blast2Go cf. enrichement.R)
{
## Transcription factors, kinases, phosphatases, chromatin remodelers
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
load("regulators.Rdata", .GlobalEnv)
lsos(pat="GO2gene_.*|.*_matrix|.*2fold|TF16")

#GOx3 <- data.frame(V1=TF2fold$locus)		## 61 TF
GOx3 <- data.frame(V1=TF16[,2])			## 17 TF
cat("\nDimentions of GOx3 rbind 3 GO sets before",dim(GOx3),"\n")
}
TFGO=GOx3
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
genes_SIG <- genes_coc[genes_coc$locus %in% TFGO$V1,]
genes_mat <- as.matrix(genes_SIG[,3:17])

## Transform & Normalize
Tgo <- normalize.quantiles(genes_mat)		# IMP: matrix where columns are arrays, rows are probes

rownames(Tgo) <- genes_SIG$locus
		#rownames(Tgo) <- paste(genes_SIG$locus,"|",genes_SIG$description)
Tgo <- t(Tgo)
cat("Dimentions of the Transcription Factor Matrix :",dim(Tgo),"\n")	#15x1444

speciesF <- as.longitudinal(Tgo, repeats=3, time=c(0,6,12,22,32))
TF <- dataFormat(speciesF)		# Transcription factors
}

==================================RUN(a)
## ANOVA ONE WAY RESULTS
## LOAD 3 subsets (from 5 different Templates)
{
setwd("C:/Dropbox/Workshop2013/Work/R/ANOVA/")
load("anova.Rdata", .GlobalEnv)
lsos(pat="nova")
}
## choose subset
N=setAnova
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T, row.names=3)
description <- genes[rownames(genes) %in% colnames(N),]
#description <- description[ ! rownames(description) %in% TFGO$V1,]
#cat("Dimentions of the Description dataframe :",dim(description),"\n")	# mx15

NNmatrix <- as.matrix(t(N))
#cat("Dimentions of the Observed Genes Before GO :",dim(NNmatrix),"\n")	# mx15

#NNsig <- NNmatrix[! rownames(NNmatrix) %in% TFGO$V1,]
cat("Dimentions of the Observed Genes :",dim(NNmatrix),"\n")	# mx15

## Transform & Normalize
NNnormalize <- normalize.quantiles(NNmatrix) 		# IMP: matrix where columns are arrays, rows are probes
	# rownames(NNnormalize) <- rownames(NNsig)
rownames(NNnormalize) <- rownames(description)
	# rownames(NNnormalize) <- paste(rownames(description),"|",description$description, sep="")
NNnormalize <- t(NNnormalize)
species <- as.longitudinal(NNnormalize, repeats=3, time=c(0,6,12,22,32))
G <- dataFormat(species)		# Observed Genes
}
==================================RUN(b)
## SUBSET E
## (GROWTH, DEVELOPMENT, LOCOMOTION, RHYTHMIC)
## LOAD SETE
{
## SUBSET E GROWTH, DEVELOPMENT, LOCOMOTION, RHYTHMIC
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
load("subsetE2.Rdata", .GlobalEnv)
lsos(pat="GO2gene_.*|.*_matrix|.*setE")
## view enrichment.R for more detail

setEx3 <- data.frame(V1=subsetE_2fold2[,2])	## neural, tissue, muscle
#setEx3 <- data.frame(V1=subsetE_2fold2[,2])	## growth, development, locomotion, rhythmic
cat("\nDimentions of setEx3 rbind 3 Processes sets before",dim(setEx3),"\n")
}
x3setE=setEx3
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansW2.txt", sep="\t", header=T, row.names=3)
description <- genes[rownames(genes) %in% x3setE$V1,]
cat("Dimentions of the Observed Genes before TF extraction :",dim(description),"\n")	# mx15
description <- description[ ! rownames(description) %in% TFGO$V1,]
description2 <- as.matrix(description[,-c(1,2)])
## Transform & Normalize
NNnormalize <- normalize.quantiles(description2) 		# IMP: matrix where columns are arrays, rows are probes
cat("Dimentions of the Observed Genes After TF extraction :",dim(NNnormalize),"\n")	# mx15
	# rownames(NNnormalize) <- rownames(NNsig)
rownames(NNnormalize) <-rownames(description)
NNnormalize <- t(NNnormalize)
species <- as.longitudinal(NNnormalize, repeats=3, time=c(0,6,12,22,32))
G <- dataFormat(species)		# Observed Genes
}

==================================RUN(c)
## CLUSTERS (12 clusters)
## Load all clusters 12 from (last soft clustering for paper 2 w all.set -- Library 1) (optional)
{					## Load SoftClustering Mfuzz package --- 12 clusters
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
load("clusters12_Lm.Rdata", .GlobalEnv)
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
load("enrichment.Rdata", .GlobalEnv)
lsos(pat="cluster.*")
}
## Clusters (should load unique GO first)		------ choose either N -------
N=cluster1
N=data.frame(x=cluster12_fisher_alpha07_sig)		## enriched DETs
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
N <- genes_coc[genes_coc$locus %in% N$x,]
cat("Before cleaning :", dim(N),"\n")

## EXTRACT genes WITH 2 FOLD CHANGE OR MORE
N2 <- NULL
for(i in 3:ncol(N)){
N2 <- rbind(N2, N[N[,i] >=1,])
}
N2 <- N2[!duplicated(N2$locus),]
cat("Two fold at any stage :", dim(N2),"\n")
## Remove TF (transcription factors) from all.set
Nsig <- N2[! N2[,2] %in% TFGO$V1,]
cat("After TF cleaning :", dim(Nsig),"\n")

#N <- as.matrix(N[,3:17])
N <- as.matrix(Nsig[,3:17])
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
==================================RUN(c)
## SEARCH
## Search and locate where TF (extracted from Library 1) in GOx3 are in clusters of Library 1 (12 clusters)
{						## RUN (1)
clusters <- list(cluster1,cluster2,cluster3,cluster4,cluster5,cluster6,cluster7,cluster8,cluster9,cluster10,cluster11,cluster12)
Lo=as.vector(as.matrix(GOx3))
Lo
}
## SET (LoX)
LoX=c("Locus29915")
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




=======
## START
=======
tmp <- runif(1); set.seed(426924)
## Calculate Hankel Hidden States at Lag=1, 2, 3
K <- hankel(G,lag=1,cutoff=0.9,type="median")$dim ; K

## Feedback networks with/without hidden states
networkB1 <- ebdbn(input="Feedback",G,K=3,conv.1=0.15, conv.2=0.01, conv.3=0.001, verbose=T)
networkD1 <- ebdbn(input="Feedback",TF,K=3,conv.1=0.15, conv.2=0.01, conv.3=0.001, verbose=T)
gc()
## Input networks with/without hidden state
network <- ebdbn(input=TF,G,K=3, conv.1=0.15, conv.2=0.05, conv.3=0.001, verbose=T)
networkH1 <- ebdbn(input=TF,G,K=4, verbose=T)
networkG12 <- ebdbn(input=TF,G,K=3, conv.1=0.15, conv.2=0.01, conv.3=0.001, verbose=T)
plot(networkG12, sig.level=0.99, interactive=F, clarify=T,layout = layout.fruchterman.reingold)
gc()
## END


==================================
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

===========================
## SAVED NETWORKS AND DESCRIPTIONS
===========================
## WORKING (BUT USING OLD MA PREMPROCESSING, with normexp, w betweenNormalization)
*net6 : feedback | EGG | vs GOx3 | 0001 | k=3 | .99
*net7 : feedback | VELI | vs GOx3 | 01 005 0001 | k=3 | .95
*net21 : feedback | setE | vs GOx3 | 0001 | k=1 | .99
**net14 : input | EGG | vs GOx3 | 015 005 001 | k=3 | .99
**net15 : input | VELI | vs GOx3 | 015 001 001 | k=3 | .99
**net16 : input | JUV | vs GOx3 | 01 005 0001 | k=3 | .99
**net20 : input | setE | vs GOx3 | 01 005 0001 | k=3 | .99
*net30 : feedback | setE | 01 005 0001 | k=3 | .98
## WORKING NEW MA PREPROCESSING W BG CORRECT, WEIGHTED, AND W/O BETWEEN NORMALIZATION
*networkB : feedback | all anova | 015 001 0001 | k=3 | .99
networkC : feedback | sig_All_matrix | 015 001 0001 | k=3 | .99 		## TF that are found in Library 1 (186)
*networkD : feedback | TF2fold | 015 001 0001 | k=3 | .99			## TF found in library 1 and have a minimum of 2 fold change  (61)
*networkE : input | TF2fold | vs All anova | 015 001 0001 | k=3 | .99
**networkE1 : input | TF2fold | vs All anova | 015 001 0001 | k=3 | .99
networkF : input | subsetE_2fold | TF17 | 015 001 0001 | k=3 | .99
# Final new set
****networkB1 : feedback | all anova | 015 001 0001 | k=3 | .99
****networkD1 : feedback | TF2fold | 015 001 0001 | k=3 | .99			## TF found in library 1 and have a minimum of 2 fold change  (61)
****networkE3 : input | TF2fold | vs All anova | 015 005 0001 | k=3 | .99
****networkE4 : input | TF16 | vs All anova | 015 005 0001 | k=3 | .99 |
****networkF2 : input | subsetE_2fold | TF2fold | 01 001 0001 |  k=4 | .99
****networkF1 : input | subsetE_2fold | TF16 |  k=4 | .99		## dont use Gs and TFs in cytoscape they are false (only w/ TF16)
****networkG : input | clusters (enriched w MF) | TF16 |  015 001 0001 | k=3 | .99		## dont use Gs and TFs in cytoscape they are false (only w/ TF16)
****networkH1 : input | subsetE_2fold2 | TF16 |  k=4 | .99		## dont use Gs and TFs in cytoscape they are false (only w/ TF16)

=======================
## SAVE A COPY OF THE MEMORY
=======================
setwd("C:/Dropbox/Workshop2013/Work/R/ebdbn/")
lsos(pat="network.*")
save(list=ls(pattern="network.*"),file="Networks2.3.Rdata")
load("Networks2.3.Rdata", .GlobalEnv)

#########################################################################
#								SECTION -2-								#
#########################################################################
## EXTRACT NETWORK INTERATION TO CYTOSCAPE

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
## CHOOSE NETWORK
NET = networkG1
plot(NET, sig.level=0.99, interactive=F, clarify=T)

## FUNCTION THAT TAKE THE ZSCORES $Z AND CREATE A BINARY MATRIX WITH 1 AND -1 FOR UP AND DOWN REGULATED DET FOR 0.99 THRESHOLD
net2Bin <- function(x){

for(i in 1:length(x))
		if(x[i]>='2.57'){
		x[i]=1
			} else{
			if(x[i]>'0'){
			x[i]=0
				} else{
				if(x[i]>='-2.57'){
				x[i]=1
					} else{
					if(x[i]>'-0'){
					x[i]=0}}
}
}
x
}
sigBinary	 <-net2Bin(NET$z)

## FUNCTION TO EXTRACT REGULATION INTO CYTOSCAPE FORMAT (FROM RAU)
visualizeNet <- function (zscores, sig, type = "feedback")  {
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

## VISUALIZE EXTRACTS CONNECTIONS --- (choose Inp or Feed)
cytoNet <- visualizeNet(NET$z, sigBinary, type="input")
cytoNet <- visualizeNet(NET$z, sigBinary, type="feedback")
## EXTRACT INTERACTIONS
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\ebdbn\\Second\\")
write.table(cytoNet, "networkG12.txt",sep="\t", quote=F, col.names=F, row.names=F)


{	## EXTRACT GENE EXPRESSION (change .txt)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
NN <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)	## W/n+Weight
require(plyr)
N <- t(NN[,-c(1,2)])
N <- data.frame(X1=sort(rep(seq(1:5),3)),N)
N <- ddply(N, .(X1), numcolwise(mean))
N <- t(N[,-1])
rownames(N) <- rownames(NN)
colnames(N) <- c("Egg","Troco","Veli","Pedi","Juv")
head(N)
N <- data.frame(locus=NN[,2],ID=NN[,1],N)
}
{	## EXPORT GENE EXPRESSION TO CYTOSCAPE TABLES
## TF1 + logs1 ... (change file name)
NlogsTF <- N[N$locus %in% colnames(NET$z),] ; dim(NlogsTF)
nbAc=ncol(NET$z)
TFlogs <- data.frame(X1=paste("TF",seq(1:nbAc),sep=""),NlogsTF)
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\ebdbn\\Second\\")
write.table(TFlogs, "networkG1_TF_logs.txt",sep="\t", quote=F, col.names=F, row.names=F)

## G1 + logs1 ... (change file name)
NlogsG <- N[N$locus %in% rownames(NET$z),] ; dim(NlogsG)
nbAr=nrow(NET$z)
Glogs <- data.frame(X1=paste("G",seq(1:nbAr),sep=""),NlogsG)
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\ebdbn\\Second\\")
write.table(Glogs, "networkG12_G_logs.txt",sep="\t", quote=F, col.names=F, row.names=F)
}
{	## EXTRACT NAMES OF TF AND GENES (optional)
## TF1, TF2 .... (change file name)
nbAc=ncol(NET$z)
TFlocus <- data.frame(X1=paste("TF",seq(1:nbAc),sep=""),X2=colnames(NET$z))
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\ebdbn\\")
write.table(TFlocus, "networkF_TF.txt",sep="\t", quote=F, col.names=F, row.names=F)

## G1, G2 .... (change file name)
nbAr=nrow(NET$z)
Glocus <- data.frame(X1=paste("G",seq(1:nbAr),sep=""),X2=rownames(NET$z))
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\ebdbn\\")
write.table(Glocus, "networkF_G.txt",sep="\t", quote=F, col.names=F, row.names=F)
}
## INPUT NETWORKS ONLY
=================
{	## TF WITH TRUE BOOLEAN APPEND TO LOCUS (TO CHANGE THE NODE TYPE IN CYTOSCAPE)
nbAc=ncol(NET$z)
TFtrue <- data.frame(X1=paste("TF",seq(1:nbAc),sep=""),X2="true")
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\ebdbn\\Second\\")
write.table(TFtrue, "networkG1_TF_true.txt",sep="\t", quote=F, col.names=F, row.names=F)
}

#########################################################################
#								SECTION -2-								#
#########################################################################
## READ FILE CONTAINING TF FITTED TO ANOVA SET. BEST TF HIGH LOGS AND HIGH CONNECTIVITY (SIGNIFICANT HUBS)
TF17 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ebdbn\\networkE_17TF.txt")
lsos(pat="GO2gene_.*|.*_matrix|.*2fold|TF17")


#########################################################################
#								SECTION -3-								#
#########################################################################
## CYTOSCAPE VISUALS (MOST INTERESTING)

Run "Network Analyser"

edge line type / interaction / discrete
edge stroke color (unselected) / EdgeBetwenness / continuous
node fill colors / Egg / continuous
node size / averageShortesPAth / Continuous

node fill colors / AverageShortPath / discrete / rainbow
node label font size / AverageShortPath / discrete / number series / 3 / 3
edge width/ edgeBetweeness / continuous
edge stroke color (unselected) / interaction / discrete



#########################################################################
#								SECTION -4-								#
#########################################################################

library(ebdbNet)
tmp <- runif(1) ## Initialize random number generator
set.seed(935488) ## Set seed
{
P <- 10 ## 10 genes
## Create artificial true D matrix
{
Dtrue <- matrix(0, nrow = P, ncol = P)
index <- expand.grid(seq(1:P),seq(1:P))
selected.index <- sample(seq(1:(P*P)), ceiling(0.25 * P * P))
selected.edges <- index[selected.index,]
for(edge in 1:ceiling(0.25 * P * P)) {
tmp <- runif(1)
if(tmp > 0.5) {
Dtrue[selected.edges[edge,1], selected.edges[edge,2]] <-runif(1, min = 0.2, max = 1)
}
else {
Dtrue[selected.edges[edge,1], selected.edges[edge,2]] <-runif(1, min = -1, max = -0.2)
}
}
}
## Create artificial estimated D matrix
{
Dest <- matrix(0, nrow = P, ncol = P)
index <- expand.grid(seq(1:P),seq(1:P))
selected.index <- sample(seq(1:(P*P)), ceiling(0.25 * P * P))
selected.edges <- index[selected.index,]
for(edge in 1:ceiling(0.25 * P * P)) {
tmp <- runif(1)
if(tmp > 0.5) {
Dest[selected.edges[edge,1], selected.edges[edge,2]] <-runif(1, min = 0.2, max = 1)
}
else {
Dest[selected.edges[edge,1], selected.edges[edge,2]] <-runif(1, min = -1, max = -0.2)
}
}
}
SS <- calcSensSpec(Dtrue, Dest)
}

SS <- calcSensSpec(NET$z, NET$DPost)
(SS$TP/(SS$TP+SS$FN))
(SS$TN/(SS$TN+SS$TP))

tmp <- runif(1) ## Initialize random number generator
## Generate artificial values for sensitivity and complementary specificity.
fn <- function(x) {return(-1/(x^7)+1)}
set.seed(6)
sens <- c(fn(seq(1, 2.7, length = 100)),1) ## Sensitivity
cspec <- seq(0, 1, by = 0.01) ## Complementary specificity
## Calculate the AUC of the ROC curve
 calcAUC(sens, cspec) ## AUC of this ROC curve is 0.9030868

sens <- c(0.2, 0.24,0.28)
cspec <- c(0.9167, 0.9032,0.8906)

#########################################################################
#								SECTION -5-								#
#########################################################################
## EXTRACT TABLES
library(xtable)
helpWlatex <- read.table("C:\\Users\\Sleiman Bassim\\Downloads\\AllnetworkG.txt", sep="\t",header=T)
head(helpWlatex)
table2latex <- xtable(helpWlatex)
print(table2latex, floating=T)


#########################################################################
#								SECTION -6-								#
#########################################################################
In calculating ROC curves for networks, you have to have 2 things:

1) The true network -- in the case of the ebdbNet documentation for calcSensSpec, that corresponds to Dtrue. For the purposes of illustrating how to use the calcSensSpec function, I took Dtrue to be just an arbitrary adjacency matrix (which of course doesnt really make sense in practice). What you would use for Dtrue in practice is the adjacency matrix that corresponds to the true underlying graph (with 1's in positions corresponding to true edges, and 0's otherwise), for example the graph structure you use to simulate data. If you are using real data without a known graph, there would of course be no way to calculate ROC curves (how would you know what are true positives, false positives, etc.?!).

2) The estimated network -- this corresponds to the adjacency matrix of the estimated graph (with 1's in positions corresponding to true edges, and 0's otherwise) at a range of different significance thresholds. So if the significance threshold is very small, you'll have almost all 1's whereas if it is large, you'll have almost all 0's. Again, in the calcSensSpec function documentation, in the interest of illustrating how to call the function, I took Dest to be just an arbitrary adjacency matrix (again, nonsensical in practice).

Below you can find a toy example showing how to calculate the AUROC (and sensitivity and specificity over a range of thresholds) from the ebdbNet results; that may hopefully help clarify things a bit. Finally, to do a similar (and simplified!) calculation using the pROC package, you would use the same 0-1 true data matrix, and the absolute value of the z-scores from the estimated matrix -- an example is included at the bottom of the code below.

Hope that helps!
Best, Andrea

library(ebdbNet)
library(pROC)
tmp <- runif(1) ## Initialize random number generator
set.seed(125214) ## Save seed

## Simulate data
R <- 5
T <- 10
P <- 10
simData <- simulateVAR(R, T, P, v = rep(10, P), perc = 0.10)
Dtrue <- simData$Dtrue
y <- simData$y

####################################################
## Run EB-DBN without hidden states
####################################################
## Choose alternative value of K using hankel if hidden states are to be estimated
K <- hankel(y)$dim
## Run algorithm
net <- ebdbn(y = y, K = 0, conv.1 = 0.15, conv.2 = 0.10, conv.3 = 0.10,
    verbose = TRUE)
plot(net, sig.level=0.99, interactive=F, clarify=T,layout = layout.fruchterman.reingold)
####################################################
## Calculate sensitivity and specificity at various significance thresholds
## Then calculate AUROC
####################################################
sens <- spec <- rep(NA, length = length(seq(0,1,by=0.01)))
index <- 1
for(alpha in seq(0,1,by=0.01)) {
    cutoff <- qnorm((1 + alpha)/2)
    Dest <- matrix(0, nrow=nrow(net$z), ncol=ncol(net$z))
    Dest[which(abs(net$z) > cutoff, arr.ind = TRUE)] <- 1
    tmp <- calcSensSpec(trueMatrix=Dtrue, estMatrix=Dest)
    sens[index] <- tmp$TP / (tmp$TP + tmp$FN)
    spec[index] <- tmp$TN / (tmp$TN + tmp$FP)
    index <- index + 1
}
AUC <- calcAUC(sens, 1-spec)
## Plot ROC curve
plot(1-spec, sens, type="l")

####################################################
## Calculate AUROC using pROC package
####################################################
roc.calc <- roc(Dtrue, abs(net$z))
roc.calc$auc ## Not quite the same, different approximation used
## Plot ROC curve
plot(roc.calc)


roc.calc <- roc(sigBinary, abs(NET$z))





