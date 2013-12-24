
#########################################################################
#																			SECTION -1-													#
#########################################################################
# ENRICHMENT ANALYSIS FROM( ALEXA2010TOPGO)
# FIRST STEP - CREATE A GENE-TO-GO OBJECT (BY EXTRACTING GO-TO-GENE FROM BLAST2GO)
# SECOND STEP - PREDEFINE A LIST OF INTERESTING GENES (LIBRARY 1)

{					## Load L.Rdata
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
load("enrichment.Rdata", .GlobalEnv)
lsos(pat="")
}

library(topGO)
## choose from BP ; MF ; CC
x="GO2Locus_mf.txt"
y=2							## cluster 1-12 ALL
z=cluster11$x		## cluster 1-12 but only >0.7

## Build Custom GO annotation File from the 15000 annotated transcripts (Blast2GO annotations)
{
## x is the file extracted from Blast2GO where every GOterm is associated to available genes
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
GO2geneID <- readMappings(x)
str(head(GO2geneID))
## Every Locus with its GO term
geneID2GO <- inverseList(GO2geneID)
str(head(geneID2GO))
length(geneID2GO)		## number of transripts annotated
}
## Extract Library 1 gene-to-GO
## OPTION 1 (all library 1)
## allSet_lm is Library 1
{
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
all.set <- read.table("allSet_Lm.txt", sep="\t", header=T)
str(head(all.set))
dim(all.set)
}
## OPTION 2 (only high membership genes of clusters of Library 1)
## Load all clusters from (last soft clustering for paper 2 w all.set)
{					## Load SoftClustering Mfuzz package --- 12 clusters
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
load("clusters12_Lm.Rdata", .GlobalEnv)
lsos(pat="cluster.*")
}
## OPTION 3 (All genes of the clusters of Library 1)
lsos(pat="acor.*")
## OPTION 4 , ANOVA RESULTS
{
setwd("C:/Dropbox/Workshop2013/Work/R/ANOVA/")
load("anova.Rdata", .GlobalEnv)
lsos(pat="nova")
}

=======
## START
=======
## choose which genes to consider
{			## change BP - MF -C
geneNames <- names(geneID2GO)
head(geneNames)
length(geneNames)
## choose :
	#geneList <- factor(as.integer(geneNames %in% all.set[,1]))		## Library 1
	#geneList <- factor(as.integer(geneNames %in% colnames(setAnova)))		## Anova
	#geneList <- factor(as.integer(geneNames %in% acore_list01[[y]][,1]))		## all clustered genes
	geneList <- factor(as.integer(geneNames %in% z))		## only genes w high membership values
names(geneList) <- geneNames
str(geneList)
## Build TopGOdata object (before enrichment analysis
GOdata <- new("topGOdata", description = "Library 1", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO=geneID2GO)
GOdata
}
## compute tests
{
resultFis <- runTest(GOdata, algorithm="classic", statistic="fisher")
resultKS <- runTest(GOdata, algorithm="weight01", statistic="ks")
result.elim <- runTest(GOdata, algorithm="elim", statistic="ks")
}
## Extract Results
cluster12_fisher_alpha07 <- GenTable(GOdata, classicFisher=resultFis,orderBy="classicFisher",topNodes=14)
cluster12_fisher_alpha07
allRes <- GenTable(GOdata, classicFisher=resultFis, KS=resultKS, elimKS=result.elim ,orderBy="classicFisher",ranksOf="classicFisher",topNodes=20)
allRes

===============
## ADDITIONAL TOOLS
===============
## SEARCH FOR LOCUS WITH SPECIFIC GO TERM
{		goID="GO:0044444"
goID_df <- as.matrix(GO2geneID[goID][[1]])
sort(goID_df[goID_df %in% colnames(setAnova),])
}
## SEARCH FOR GO term WITH SPECIFIC LOCUS
{		locus="Locus33502"
 as.matrix(geneID2GO[locus][[1]])
}
## SEARCH FOR ALL LOCUS OF A SET
locus_df <- as.matrix(geneID2GO[names(setAnova)])

## OPTIONAL
{
## Enrichment analysis
whichTests()
whichAlgorithms()
## p-value extraction
{
pvalFis <- score(resultFis)
pval.elim <- score(result.elim)[names(pvalFis)]
pvalKS <- score(resultKS)[names(pvalFis)]
}
## Show significant
geneData(resultFis)
geneData(resultKS)
geneData(result.elim)
## Histograms
hist(pvalFis, 50, xlab="p-values fisher")
hist(pval.elim,50,xlab="p-values ks/elim")
hist(pvalKS, 50,xlab="p-values ks")
## plot
{
gstat <- termStat(GOdata, names(pvalFis))
gSize <- gstat$Annotated / max(gstat$Annotated)*4
plot(pvalFis, pval.elim, xlab="p-value classic", ylab="p-value ks/elim", pch=19, cex=gSize)
}
}

## Obtain the list of genes
## significant genes by GO-terms
{
cluster12_fisher_alpha07_sig <- sigGenes(GOdata)
str(cluster12_fisher_alpha07_sig)
termStat(GOdata)
}

## Save a copy of the memory
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
lsos(pat="cluster.*_")
save(list=ls(pattern="cluster.*_"),file="enrichment.Rdata")
#info:
clusterN_fisher : all genes of a cluster (acore_list01) (MF enriched)
clusterN_fisher_alpha07 : table of cluster genes >0.7 (second round of clustering) (MF enriched)
clusterN_fisher_alpha07_sig : only genes (MF enriched)

## EXTRACT TABLES
library(xtable)
cluster_table <- xtable(cluster1_fisher)
print(cluster_table, floating=T)

## Just for testing
## Select a GO term from the indexed ones
x="GO:0005198"
{
goID <- x
## All genes
gene.universe <- genes(GOdata)
go.genes <- genesInTerm(GOdata, goID)[[1]]
## Significant genes
sig.genes <- sigGenes(GOdata)

my.group <- new("classicCount", testStatistic = GOFisherTest, name = "fisher", allMembers = gene.universe, groupMembers = go.genes,sigMembers = sig.genes)
contTable(my.group)
runTest(my.group)
}


#########################################################################
#								SECTION -2-								#
#########################################################################
## SQLite (ref pages2008annotation)
## working

## Retrieving data from a web source (Uniprot)
## use only IDs not protein names
ids <- c('P13368','Q6GZX4')
{
uri <- 'http://www.uniprot.org/uniprot/?query='
idStr <- paste(ids, collapse="+or+")
formatSQ <- '&format=tab'
	colSQ <- '&columns=id,sequence'
fullUri <- paste0(uri, idStr,formatSQ)
read.delim2(fullUri, stringsAsFactors=F)
}

## Function
getUniprot <- function(query, cols)
{
## query and cols start as a character vectors
qstring <- paste(query, collapse="+or+")
cstring <- paste(cols, collapse=",")
uri <- 'http://www.uniprot.org/uniprot/?query='
fullUri <- paste0(uri,qstring,'&format=tab&columns=',cstring)
dat <- read.delim(fullUri, stringsAsFactors=FALSE)
## now remove things that were not in the specific original query...
dat <- dat[dat[,1] %in% query,]
dat
 }


#########################################################################
#								SECTION -3-								#
#########################################################################
## EXTRACT TRANSCRIPTION FACTORS (TF); AND CHAPERONES (CH): KINASES-PHOSPHATASES-CHROMATIN REMODELERS
## LOAD THE BELOW PREPARED SETS
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
load("regulators.Rdata", .GlobalEnv)
lsos(pat="GO2gene_.*|.*_matrix|.*2fold|TF16")

## EXTRACT GO TERMS ASSOCIATED WITH "TRANSCRIPTION" ; "KINASE" ; "PHOSPHATASE" ; "CHROMATIN REMODELER"
## CHANGE Y AND FILE NAMES
y = "transcription"
temp <- readLines("C:\\Blast2GO\\Annotated15K\\GO2Locus_mf_raw.txt")
cat(temp[grep(y, temp)], file="GO_clean_remodeler.txt", sep="\n")

## EXTRACT LOCUS
## CHANGE X AND TF/P/K/C
require(topGO)
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
x="GO_clean_chromatin.txt"
GO2gene_C <- readMappings(x)
geneID2GO_C <- inverseList(GO2gene_C)
TF_matrix <- as.matrix(names(geneID2GO_TF))
K_matrix <- as.matrix(names(geneID2GO_K))
P_matrix <- as.matrix(names(geneID2GO_P))
C_matrix <- as.matrix(names(geneID2GO_C))

## MERGE AND REMOVE DUPLICATED
All_matrix <- rbind(TF_matrix, K_matrix, C_matrix, P_matrix)
All_matrix <- All_matrix[! duplicated(All_matrix)]

## extract significant TF , P, C, K
## Library 1 for subset A, B, C w background correction to minimum | w/o between Arrays normalization | weighted (-1,2)
{
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
all.set <- read.table("allSet_Lm.txt", sep="\t", header=T)
}
sig_All_matrix <- all.set[all.set$x %in% All_matrix,]
sig_TF_matrix <- all.set[all.set$x %in% TF_matrix,]
sig_K_matrix <- all.set[all.set$x %in% K_matrix,]
sig_P_matrix <- all.set[all.set$x %in% P_matrix,]
sig_C_matrix <- all.set[all.set$x %in% C_matrix,]


## EXTRACT GENE EXPRESSION (change .txt)
{
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
## EXPORT GENE EXPRESSION TO CYTOSCAPE TABLES
{
NlogsTF <- N[N$locus %in% sig_All_matrix,] ; dim(NlogsTF)
nbAc=ncol(NET$z)
networkA_TF_logs <- data.frame(X1=paste("TF",seq(1:nbAc),sep=""),NlogsTF)
}
## EXTRACT TF WITH 2 FOLD CHANGE OR MORE
TF2fold <- NULL
for(i in 4:8){
TF2fold <- rbind(TF2fold, networkA_TF_logs[networkA_TF_logs[,i] >=1,])
}
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\ebdbn\\")
write.table(TF2fold, "TF2fold.txt",sep="\t", quote=F, col.names=F, row.names=F)




## SAVE
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
lsos(pat="GO2gene_.*|.*_matrix|.*2fold|TF16")
save(list=ls(pattern="GO2gene_.*|.*_matrix|.*2fold|TF16"), file="regulators.Rdata")


#########################################################################
#								SECTION -4-								#
#########################################################################
## EXTRACT SUBSET E (first half)
temp <- readLines("C:\\Blast2GO\\Annotated15K\\GO2locus_bp_raw.txt")
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
## CHANGE Y TO GROWTH ; DEVELOPMENT ; LOCOMOTION ; RHYTHMIC
y="growth"
cat(temp[grep(y,temp)], file="subsetE_rhythmic.txt", sep="\n")


## EXTRACT LOCUS
## REMOVE DESCRIPTION FROM GREP FILES
## CHANGE X AND G/D/L/R
require(topGO)
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
x="subsetE_rhythmic.txt"
GO2gene_R <- readMappings(x)
geneID2GO_R <- inverseList(GO2gene_R)
G_matrix <- as.matrix(names(geneID2GO_G))
D_matrix <- as.matrix(names(geneID2GO_D))
L_matrix <- as.matrix(names(geneID2GO_L))
R_matrix <- as.matrix(names(geneID2GO_R))

## MERGE AND REMOVE DUPLICATED
All_matrix <- rbind(G_matrix, D_matrix, L_matrix, R_matrix); dim(All_matrix)
All_matrix <- All_matrix[! duplicated(All_matrix)]; length(All_matrix)

## EXTRACT SIGNIFICANT G, L, D, R FROM LIBRARY 1
## all.set for subset A, B, C w background correction to minimum | w/o between Arrays normalization | weighted (-1,2)
{
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
all.set <- read.table("allSet_Lm.txt", sep="\t", header=T)
}
sig_All_matrix <- all.set[all.set$x %in% All_matrix,]
sig_G_matrix <- all.set[all.set$x %in% G_matrix,]
sig_D_matrix <- all.set[all.set$x %in% D_matrix,]
sig_L_matrix <- all.set[all.set$x %in% L_matrix,]
sig_R_matrix <- all.set[all.set$x %in% R_matrix,]


## EXTRACT GENE EXPRESSION (change .txt)
{
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
## EXPORT GENE EXPRESSION TO CYTOSCAPE TABLES
logsSubsetE <- N[N$locus %in% sig_All_matrix,] ; dim(logsSubsetE)
nbAc=nrow(logsSubsetE)
logsSubsetE <- data.frame(X1=paste("G",seq(1:nbAc),sep=""),logsSubsetE)

## EXTRACT TF WITH 2 FOLD CHANGE OR MORE
subsetE_2fold <- NULL
for(i in 4:8){
subsetE_2fold <- rbind(subsetE_2fold, logsSubsetE[logsSubsetE[,i] >=1,])
}
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\ebdbn\\")
write.table(subsetE_2fold, "subsetE_2fold.txt",sep="\t", quote=F, col.names=F, row.names=F)


names(subsetE_2fold)
length(sig_R_matrix[sig_R_matrix %in% subsetE_2fold$locus])

## SAVE
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
lsos(pat="GO2gene_.*|.*_matrix|.*setE")
save(list=ls(pattern="GO2gene_.*|.*_matrix|.*setE"), file="subsetE.Rdata")


#########################################################################
#								SECTION -5-								#
#########################################################################
## EXTRACT SUBSET E (second half)
temp <- readLines("C:\\Blast2GO\\Annotated15K\\GO2locus_bp_raw.txt")
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
## CHANGE Y TO Muscle, Neural, tissue development
y="tissue"
cat(temp[grep(y,temp)], file="subsetE_muscle.txt", sep="\n")
cat(temp[grep(y,temp)], file="MFtissue.txt", sep="\n")
é

## EXTRACT LOCUS
## REMOVE DESCRIPTION FROM GREP FILES
## CHANGE X AND G/D/L/R
require(topGO)
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
x="subsetE_tissue.txt"
GO2gene_Ti <- readMappings(x)
geneID2GO_Ti <- inverseList(GO2gene_Ti)
Ne_matrix <- as.matrix(names(geneID2GO_Ne))	## neural
Mu_matrix <- as.matrix(names(geneID2GO_Mu))	## muscle
Ti_matrix <- as.matrix(names(geneID2GO_Ti))		## tissue

## MERGE AND REMOVE DUPLICATED
All_matrix2 <- rbind(Ne_matrix, Ti_matrix,Mu_matrix); dim(All_matrix2)
All_matrix2 <- All_matrix2[! duplicated(All_matrix2)]; length(All_matrix2)

## EXTRACT GENE EXPRESSION (change .txt)
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
NN <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)	## W/n+Weight
}
## EXPORT GENE EXPRESSION TO CYTOSCAPE TABLES
logsSubsetE2 <- NN[NN$locus %in% All_matrix2,]
cat("All of extracted entries :", dim(logsSubsetE2), "\n")

## EXTRACT TF WITH 2 FOLD CHANGE OR MORE
subsetE_2fold2 <- NULL
for(i in 3:ncol(NN)){
subsetE_2fold2 <- rbind(subsetE_2fold2, logsSubsetE2[logsSubsetE2[,i] >=1,])
}
cat("After 2fold clean :", dim(subsetE_2fold2), "\n")
subsetE_2fold2 <- subsetE_2fold2[ ! duplicated(subsetE_2fold2$locus),]
cat("After removal of redundant entries:", dim(subsetE_2fold2), "\n")

## SAVE
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
lsos(pat="GO2gene_.*|.*_matrix|.*setE")
save(list=ls(pattern="GO2gene_.*|.*_matrix|.*setE"), file="subsetE.Rdata")



