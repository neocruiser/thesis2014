#########################################################################
#								SECTION -1-								#
#########################################################################
## SVD & PCA


## Load ALL GENES
## PLOT ALL genes
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
expressed <- as.matrix(genes[,3:17])
logs <- as.matrix(t(expressed))
## Centring & Weigthed SVD
## Matrix is multiplied by a constant (standardization)
species <- NULL
for(j in 1:15) species <- rbind(species, (logs[j,]-(apply(logs,2,mean)))/(sqrt(nrow(logs)*ncol(logs))))
species.SVD <-svd(species)
## *********** COVARIANCE BIPLOT (col=principal, row=std coordinates) ****************
species.F <- sqrt(nrow(species))*species.SVD$u 
species.G <- sqrt(ncol(species))*species.SVD$v %*% diag(species.SVD$d)
plot(rbind(species.F[,1:2],species.G[,1:2]), type="n", asp=1,xlab="dim1",ylab="dim2", cex.axis=0.7, las=1)
abline(v=0, h=0, lty="dotted")
points(species.G[,1], species.G[,2], pch=20, col="grey", cex=0.5)
}

## Extract PATTERN DISCOVERY selected  eigenes
## SAM Analysis
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
load("sets.Rdata", .GlobalEnv)
lsos(pat="quanti_sigg.*|set.*")
z <- rainbow(30)		## create color palette
## SUBSET C
{
variaC <- as.integer(rownames(genes[genes$locus %in% subsetC[,1],]))
for(j in variaC) points(species.G[j,1],species.G[j,2], pch=17, col="#00E5A0", cex=0.55)
}
## SUBSET B
{
variaB1 <- as.integer(rownames(genes[genes$locus %in% setB_lo_eig1[,1],]))
for(j in variaB1) points(species.G[j,1],species.G[j,2], pch=17, col="#04B8E2", cex=0.45)
variaB2 <- as.integer(rownames(genes[genes$locus %in% setB_lo_eig2[,1],]))
for(j in variaB2) points(species.G[j,1],species.G[j,2], pch=15, col="#F30021", cex=0.4)
variaB3 <- as.integer(rownames(genes[genes$locus %in% setB_lo_eig3[,1],]))
for(j in variaB3) points(species.G[j,1],species.G[j,2], pch=3, col="#B0FA00", cex=0.45)
variaB4 <- as.integer(rownames(genes[genes$locus %in% setB_up_eig1[,1],]))
for(j in variaB4) points(species.G[j,1],species.G[j,2], pch=17, col="#04B8E2", cex=0.45)
variaB5 <- as.integer(rownames(genes[genes$locus %in% setB_up_eig2[,1],]))
for(j in variaB5) points(species.G[j,1],species.G[j,2], pch=15, col="#F30021", cex=0.4)
variaB6 <- as.integer(rownames(genes[genes$locus %in% setB_up_eig3[,1],]))
for(j in variaB6) points(species.G[j,1],species.G[j,2], pch=3, col="#B0FA00", cex=0.45)
}
## SUBSET A
{
variaA1 <- as.integer(rownames(genes[genes$locus %in% setA_lo[,1],]))
for(j in variaA1) points(species.G[j,1],species.G[j,2], pch=20, cex=0.4)
variaA2 <- as.integer(rownames(genes[genes$locus %in% setA_up[,1],]))
for(j in variaA2) points(species.G[j,1],species.G[j,2], pch=20, cex=0.4)
}
text(species.F[,1:15], labels=rep(c("E","T","V","P","J"),rep(3,5)),col="forestgreen", cex=1, font=4)

## SCREE PLOT
## eiganvalues
{
species.SVD$d^2
## sum of squares (SS) decomposed by the SVD = Total variance accounted for by 2D solution
## Total inertia
sum(species.SVD$d[1]^2,species.SVD$d[2]^2)*100		## cumulative percentage for the first two components
sum(species.SVD$d^2)*100
## Average variance
sum(species.SVD$d^2)/(nrow(logs)-1)
## Total varaince explained by the two first axes (just sum up the first 2)
eig.p <- 100*species.SVD$d / sum(species.SVD$d)
eig.p <- eig.p[seq(15,1)]
barplot(eig.p,horiz=T,cex.axis=0.7)
}


## *** experimental
## SOFT clustering :: cluster plotting
## Get The PCA names of every gene in each cluster
## Reference first to SAM.R section 1. Cluster ************LOW
## **************RUN SAM.R section 1 first
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
## Transcipts from Clusters
genes_cl1_lo <- genes_coc[genes_coc$locus %in% cluster1_lo,]
rn1 <- genes_cl1_lo$X
for(j in rn1) points(species.G[j,1],species.G[j,2], pch=15, col="red", cex=0.5)
genes_cl2_lo <- genes_coc[genes_coc$locus %in% cluster2_lo,]
rn2 <- genes_cl2_lo$X
for(j in rn2) points(species.G[j,1],species.G[j,2], pch=15, col="forestgreen", cex=0.5)
genes_cl3_lo <- genes_coc[genes_coc$locus %in% cluster3_lo,]
rn3 <- genes_cl3_lo$X
for(j in rn3) points(species.G[j,1],species.G[j,2], pch=15, col="blue", cex=0.5)
genes_cl4_lo <- genes_coc[genes_coc$locus %in% cluster4_lo,]
rn4 <- genes_cl4_lo$X
for(j in rn4) points(species.G[j,1],species.G[j,2], pch=15, col="yellow", cex=0.5)
genes_cl5_lo <- genes_coc[genes_coc$locus %in% cluster5_lo,]
rn5 <- genes_cl5_lo$X
for(j in rn5) points(species.G[j,1],species.G[j,2], pch=15, col="maroon", cex=0.5)
genes_cl6_lo <- genes_coc[genes_coc$locus %in% cluster6_lo,]
rn6 <- genes_cl6_lo$X
for(j in rn6) points(species.G[j,1],species.G[j,2], pch=15, col="black", cex=0.5)
genes_cl7_lo <- genes_coc[genes_coc$locus %in% cluster7_lo,]
rn7 <- genes_cl7_lo$X
for(j in rn7) points(species.G[j,1],species.G[j,2], pch=15, col="violetred2", cex=0.5)
genes_cl8_lo <- genes_coc[genes_coc$locus %in% cluster8_lo,]
rn8 <- genes_cl8_lo$X
for(j in rn8) points(species.G[j,1],species.G[j,2], pch=15, col="springgreen", cex=0.5)
genes_cl9_lo <- genes_coc[genes_coc$locus %in% cluster9_lo,]
rn9 <- genes_cl9_lo$X
for(j in rn9) points(species.G[j,1],species.G[j,2], pch=15, col="wheat1", cex=0.5)
}

## Get The PCA names of every gene in each cluster
## Reference first to SAM.R section 1. Cluster *************UP
## **************RUN SAM.R section 1 first
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
## Transcipts from Clusters
genes_cl1_up <- genes_coc[genes_coc$locus %in% cluster1_up,]
rn1 <- genes_cl1_up$X
for(j in rn1) points(species.G[j,1],species.G[j,2], pch=15, col="red", cex=0.5)
genes_cl2_up <- genes_coc[genes_coc$locus %in% cluster2_up,]
rn2 <- genes_cl2_up$X
for(j in rn2) points(species.G[j,1],species.G[j,2], pch=15, col="forestgreen", cex=0.5)
genes_cl3_up <- genes_coc[genes_coc$locus %in% cluster3_up,]
rn3 <- genes_cl3_up$X
for(j in rn3) points(species.G[j,1],species.G[j,2], pch=15, col="blue", cex=0.5)
genes_cl4_up <- genes_coc[genes_coc$locus %in% cluster4_up,]
rn4 <- genes_cl4_up$X
for(j in rn4) points(species.G[j,1],species.G[j,2], pch=15, col="yellow", cex=0.5)
genes_cl5_up <- genes_coc[genes_coc$locus %in% cluster5_up,]
rn5 <- genes_cl5_up$X
for(j in rn5) points(species.G[j,1],species.G[j,2], pch=15, col="maroon", cex=0.5)
genes_cl6_up <- genes_coc[genes_coc$locus %in% cluster6_up,]
rn6 <- genes_cl6_up$X
for(j in rn6) points(species.G[j,1],species.G[j,2], pch=15, col="black", cex=0.5)
}

#points(species.G[14000,1],species.G[14000,2], pch=21, cex=0.7, bg="red")
#arrows(0,0,species.G[,1], species.G[,2],col="chocolate",length=0.1,angle=10)
#text(1.2*species.G[,1], 1.2*species.G[,2],labels=colnames(logs),col="chocolate",font=4,cex=0.5)
#segments(species.G[,1:2])

## ******** FORM BIPLOT ***************
## Biplot (row=principal, col=std coordinates)
## Calculate coordinates
{
species.F <- sqrt(nrow(species))*species.SVD$u %*% diag(species.SVD$d)
species.G <- sqrt(ncol(species))*species.SVD$v
		## Contribution biplot - rescaling of the standard coordinates -
		#species.G <- ((ncol(species))^0.25)*species.SVD$v
plot(rbind(species.F[,1:2],species.G[,1:2]), type="n", asp=1,xlab="dim1",ylab="dim2", cex.axis=0.7)	
text(species.F[,1:15], labels=rep(c("E","T","V","P","J"),rep(3,5)),col="forestgreen", cex=0.8)
points(species.G[,1], species.G[,2], pch=20, col="grey", cex=0.6)
}
#arrows(0,0,species.G[,1], species.G[,2],col="chocolate",length=0.1,angle=10)
#text(1.2*species.G[,1], 1.2*species.G[,2],labels=colnames(logs),col="chocolate",font=4,cex=0.5)



## Mark selected genes on the previous biplot
## ANOVA MEV GENES
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
anova_SIG <- read.table("anovaW2.txt", sep="\t", header=F)
genes_SIG <- genes_coc[genes_coc$locus %in% anova_SIG$V2,]
rn <- genes_SIG$X
## Transcipts from ANOVA set
for(j in rn) points(species.G[j,1],species.G[j,2], pch=15, col="orange", cex=0.4)
}

#########################################################################
#								SECTION -2-								#
#########################################################################
## PCA of Library 1
## Load data
## DATA from subset A, B, C
{

# LOAD ALL SUBSETS A B C
{
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
load("sets.Rdata", .GlobalEnv)
lsos(pat="quanti_sigg.*|set.*")
}

setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
logs <- genes[ genes[,2] %in% all.set[,1],]
rownames(logs) <- logs[,2]
logs <- t(logs[,-c(1:2)])
}

## Centring & Weigthed SVD
## Matrix is multiplied by a constant
{
species <- NULL
for(j in 1:15) species <- rbind(species, (logs[j,]-(apply(logs,2,mean)))/(sqrt(nrow(logs)*ncol(logs))))
species.SVD <-svd(species)
## *********** COVARIANCE BIPLOT (col=principal, row=std coordinates) ****************
species.F <- sqrt(nrow(species))*species.SVD$u 
species.G <- sqrt(ncol(species))*species.SVD$v %*% diag(species.SVD$d)
plot(rbind(species.F[,1:2],species.G[,1:2]), type="n", asp=1,xlab="dim1",ylab="dim2", cex.axis=0.7, las=1)
abline(v=0, h=0, lty="dotted")
points(species.G[,1], species.G[,2], pch=1, col="grey", cex=0.7)
text(species.F[,1:15], labels=rep(c("E","T","V","P","J"),rep(3,5)),col="forestgreen", cex=0.8, font=4)
}


## SCREE PLOT
## eiganvalues
{
species.SVD$d^2
## sum of squares (SS) decomposed by the SVD = Total variance accounted for by 2D solution
## Total inertia
sum(species.SVD$d[1]^2,species.SVD$d[2]^2)*100		## for the first two components
sum(species.SVD$d^2)*100
## Average variance
sum(species.SVD$d^2)/(nrow(logs)-1)
## Total varaince explained by the two first axes (just sum up the first 2)
eig.p <- 100*species.SVD$d / sum(species.SVD$d)
eig.p <- eig.p[seq(15,1)]
barplot(eig.p,horiz=T,cex.axis=0.7)
}





