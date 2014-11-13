

#########################################################################
#		SECTION -1-								#
#########################################################################
## SINGLE METHOD USING MADE4 PACKAGE
vignette(package="made4", "introduction")

## REFER TO VIGNETTE
library(ade4)
library(made4)

## PAPER 2 & 3
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
#genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)	# PAPER 2
genes <- read.table("meansPaper3.txt", sep="\t", header=T)	# PAPER 3
names(genes)
expressed <- as.matrix(genes[,3:dim(genes)[2]])
windows()
overview(expressed)

ex.ca <- ord(expressed, type="coa")
plot(ex.ca, classvec=colnames(expressed), arraycol=c("green","blue","red","orange","grey"), genecol="pink")


#########################################################################
#	SECTION -2-								#
#########################################################################
## CA (from greenacre)
## Can accept 0 values of fatty acids
## Work with 3 replicates COCktail
## Fatty acids dataset (********no trocophore, no sums*********)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
acids <- read.table("ag.txt", sep="\t", header=T)
fat <- as.matrix(acids[c(1:6,10:12,16:18),c(28:37,39,41:47)])			# raw data
fat_raw <- as.matrix(acids[c(1:6,10:12,16:18),c(7:16,18,20:26)])		# Percentage

## Computation of CA by hand
{
P <- fat/sum(fat)
rm <- apply(P, 1, sum)
cm <- apply(P, 2, sum)
Z <- diag(sqrt(1/rm)) %*% (as.matrix(P)-rm %*% t(cm)) %*% diag(sqrt(1/cm))
svdZ <- svd(Z)
fat_F <- diag(1/sqrt(rm)) %*% svdZ$u  %*% diag(svdZ$d)
fat_G <- diag(1/sqrt(cm)) %*% svdZ$v

plot(rbind(fat_F, fat_G), type = "n", xlab = "Dim 1", ylab = "Dim 2", asp = 1, cex.axis=0.7)
abline(h=0,v=0,lty="dotted")
text(fat_F, labels = rownames(fat), col = "forestgreen", font = 2, cex = 0.7)
text(fat_G, labels = colnames(fat), col = "chocolate", font = 2, cex = 0.8)
}
#windows(width=50,height=25)
#par(mfrow=c(1,3),cex.axis=0.9,cex=0.9)

## Computation of CA using ca package
## *******scenario 1**********
{
library(ca)
plot(ca(fat), map="rowprincipal",col=c("forestgreen", "chocolate"))
summary(ca(fat))
}

## There is a difference between raw, percentage and ranged plotted data
## ***********scenario 2**********
{
require(vegan)
fat <- decostand(fat, method="range")
plot(ca(fat), map="rowprincipal",col=c("forestgreen", "chocolate"))
}

#########################################################################
#								SECTION -3-								#
#########################################################################
## ALL GENES
## Chapter 8
## cases = col, and species = rows
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansW2.txt", sep="\t", header=T)
expressed <- as.matrix(genes[,4:18])
}

## Decomposition with CA function
{
library(vegan)
library(ca)
genes.ranged <- decostand(expressed, method="range")
genes.ca <- ca(genes.ranged)
genes.F <- genes.ca$rowcoord %*% diag(genes.ca$sv)
genes.r.cont <- genes.ca$rowmass * (genes.F[,1]^2 + genes.F[,2]^2) / sum(genes.ca$sv[1:2]^2)

genes.names <- rownames(genes)
genes.names[genes.r.cont <0.001] <- "."

genes.rowsize <- log(1+exp(1)*genes.r.cont^0.3)
genes.rowsize[genes.r.cont<0.0005] <- 1

FF <- genes.ca$rowcoord
GG <- genes.ca$colcoord %*% diag(genes.ca$sv)

plot(c(FF[,1],GG[,1]), c(FF[,2],GG[,2]), type = "n", xlab = "", ylab = "", asp = 1, cex.axis=0.7)
text(FF[,1], FF[,2], labels = genes.names, cex = genes.rowsize, col = "grey")
}
#text(GG[,1], GG[,2], labels = colnames(genes), cex = 0.8, col = "forestgreen", font = 2)

## Mark selected genes on the previous biplot
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
anova_SIG <- read.table("anovaW2.txt", sep="\t", header=F)
genes_SIG <- genes_coc[genes_coc$locus %in% anova_SIG$V2,]
rn <- genes_SIG$X
## Transcipts from ANOVA set
for(j in rn) points(FF[j,1],FF[j,2], pch=20, bg="orange", cex=0.3)
}


## Contribution biplot
genes.ca$rownames <- genes.names
plot(genes.ca, map = "colgreen", mass = c(1,0), col = c("chocolate","forestgreen"), pch = c(20,24,16,1), cex=0.)
segments(0, 0, genes.ca$rowcoord[,1]*sqrt(genes.ca$rowmass), genes.ca$rowcoord[,2]*sqrt(genes.ca$rowmass), col = "yellow", cex=0.7)


#########################################################################
#								SECTION -4-								#
#########################################################################
## Reduced dimnsion manipulation (.ca$rowcoord = _F[,1] ; .ca$sv = svd$d ; .ca$colcoord = _G[,1])
## Work with 3 replicates Cocktail microarrays
## Genes
## Extract significant locus (anova) from original dataset meansW2
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
anova_SIG <- read.table("anovaW2.txt", sep="\t", header=F)
genes_SIG <- genes_coc[genes_coc$locus %in% anova_SIG$V2,]
logs <- as.matrix(t(genes_SIG[,4:18]))
dim(logs)	#15 x 126
}

## keep species in c for standardization (normal data set where r=cases, c=species)
{
library(ca)
library(vegan)
logs.range <- decostand(logs, method="range")

## symmetric biplot
plot(ca(t(logs.range)))
## Put species in r for ca
genes.ca <- ca(t(logs.range))
genes.ca$rowcoord[,2] <- -genes.ca$rowcoord[,2]
genes.ca$colcoord[,2] <- -genes.ca$colcoord[,2]

par(cex.axis = 0.7)
plot(genes.ca, labels = 0, col = c("forestgreen", "chocolate"))
genes.F <- genes.ca$rowcoord %*% diag(genes.ca$sv)
genes.G <- genes.ca$colcoord %*% diag(genes.ca$sv)
text(genes.F, labels = genes.ca$rownames, pos = 4, offset = 0.3, col = "forestgreen", cex = 0.8, font = 2)
text(genes.G, labels = genes.ca$colnames, pos = 4, offset = 0.3, col = "chocolate", cex = 0.8, font = 4)
text(max(genes.G[,1]), 0, "0.0571 (82.1%)", adj = c(0.6,-0.6), cex = 0.7)
text(0, max(genes.G[,2]), "0.0030 (4.4%)", adj = c(-0.1,-3), cex = 0.7)
}

## showing categories connected to their cases
{
genes.249 <- as.numeric(logs.range[,"249"])
genes.249 <- genes.249 / sum(genes.249/15)
genes.249.coord <- NULL
## From 5 to 1 with an increment of 3 (for every stage) in the loop-equation
for(j in 5:1) {
  genes.249.coord <- rbind(apply(genes.249[(3*(j-1)+1):(3*j)] * genes.ca$colcoord[(3*(j-1)+1):(3*j),], 2, sum), genes.249.coord)
  }
genes.13753 <- as.numeric(logs.range[,"13753"])
genes.13753 <- genes.13753 / sum(genes.13753/15)
genes.13753.coord <- NULL
for(j in 5:1){
  genes.13753.coord <- rbind(apply(genes.13753[(3*(j-1)+1):(3*j)] * genes.ca$colcoord[(3*(j-1)+1):(3*j),], 2, sum), genes.13753.coord)
  }


plot(genes.ca, labels = 0, col = c("forestgreen", "chocolate"))
genes.F <- genes.ca$rowcoord %*% diag(genes.ca$sv)
genes.G <- genes.ca$colcoord %*% diag(genes.ca$sv)
text(genes.F, labels = genes.ca$rownames, pos = 4, offset = 0.3, col = "forestgreen", cex = 0.8, font = 2)
text(genes.G, labels = genes.ca$colnames, pos = 4, offset = 0.3, col = "chocolate", cex = 0.8, font = 4)
text(max(genes.G[,1]), 0, "0.0571 (82.1%)", adj = c(0.6,-0.6), cex = 0.7)
text(0, max(genes.G[,2]), "0.0030 (4.4%)", adj = c(-0.1,-3), cex = 0.7)

points(rbind(genes.249.coord/15, genes.13753.coord/15), pch = 22, col = "blue")
text(rbind(genes.249.coord/15, genes.13753.coord/15), labels=c("Egg", "Troco", "Veliger", "Pedi", "Juv"), pos = 4, offset = 0.3, col = "blue", cex = 0.7, font = 4)
rownames(genes.F) <- genes.ca$rownames
segments(genes.F["249",1], genes.F["249",2], genes.249.coord[,1]/15, genes.249.coord[,2]/15, col = "violetred2")
segments(genes.F["13753",1], genes.F["13753",2], genes.13753.coord[,1]/15, genes.13753.coord[,2]/15, col="red")

#points(genes.F["249",1], genes.F["249",2], pch=24)
#segments(genes.F["249",1], genes.F["249",2], genes.G[,1], genes.G[,2])
}

## Asymmetric biplot
plot(genes.ca, map = "rowprincipal", labels = c(0,2), col = c("forestgreen", "chocolate"))
