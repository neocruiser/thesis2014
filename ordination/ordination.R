### SCRIPT 1

# plot; subset; for-loop; lm; summary(lm); rbind; apply; %in%; ddply; windows; par; 

## Litterature : Biplots in practice & Leps 2003
# additional scripts in datasets.R


## Greenakre

#########################################################################
#								SECTION -1-								#
#########################################################################
=======================
## Chapter 1 (biplots)
=======================

## Read data GENES
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
anova <- read.table("anovaW2.txt", sep="\t", header=F)
logs <- data.frame(gene=as.character(anova[,1]),
					locus=as.character(anova[,2]),
						Egg=anova[,3],
						Trocophore=anova[,5],
						Veliger=anova[,7],
						Pediveliger=anova[,9],
						Juvenile=anova[,11])
row.names(logs) <- logs$locus
logs <- logs[,-c(1,2)]
}

windows(width=50,height=25)
par(mfrow=c(2,5),cex.axis=0.9,cex=0.9)

## plot
{
plot(logs[,c(1,2)], type="n", xlab="Egg", ylab="Trocophore",las=1)
points(logs[,c(1,2)], type="p", pch=21, bg="grey")	
plot(logs[,c(1,3)], type="n", xlab="Egg", ylab="Veliger",las=1)
points(logs[,c(1,3)], type="p", pch=21, bg="grey")	
plot(logs[,c(1,4)], type="n", xlab="Egg", ylab="Pediveliger",las=1)
points(logs[,c(1,4)], type="p", pch=21, bg="grey")	
plot(logs[,c(1,5)], type="n", xlab="Egg", ylab="Juvenile",las=1)
points(logs[,c(1,5)], type="p", pch=21, bg="grey")	
plot(logs[,c(2,3)], type="n", xlab="Trocophore", ylab="Veliger",las=1)
points(logs[,c(2,3)], type="p", pch=21, bg="red")	
plot(logs[,c(2,4)], type="n", xlab="Trocophore", ylab="Pediveliger",las=1)
points(logs[,c(2,4)], type="p", pch=21, bg="red")	
plot(logs[,c(2,5)], type="n", xlab="Trocophore", ylab="Juvenile",las=1)
points(logs[,c(2,5)], type="p", pch=21, bg="red")	
plot(logs[,c(3,4)], type="n", xlab="Veliger", ylab="Pediveliger",las=1)
points(logs[,c(3,4)], type="p", pch=21, bg="yellow")	
plot(logs[,c(3,5)], type="n", xlab="Veliger", ylab="Juvenile",las=1)
points(logs[,c(3,5)], type="p", pch=21, bg="yellow")		
plot(logs[,c(4,5)], type="n", xlab="Pediveliger", ylab="Juvenile",las=1)
points(logs[,c(4,5)], type="p", pch=21, bg="green")		
					text(logs[,c(2,4)], labels=rownames(logs))
					# if i want to subset special locus
					points(subset(logs, Egg<='0' & Egg >='-1' | Juvenile<='0.7' & Juvenile>='0.4',
					select=c(1,5)), type="p", pch=21, bg="red")
}

#########################################################################
#								SECTION -2-								#
#########################################################################					
======================
##Chapter 2 (biplots) / Regressions
======================
## Read data
#preapre data for linear regression "Mev anova genes" + "Growth/Performance" vs Stages
{
#Genes
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
anova <- read.table("anovaW2.txt", sep="\t", header=F)
logs <- data.frame(gene=as.character(anova[,1]),
					locus=as.character(anova[,2]),
						Egg=anova[,3],
						Trocophore=anova[,5],
						Veliger=anova[,7],
						Pediveliger=anova[,9],
						Juvenile=anova[,11])
row.names(logs) <- logs$locus
locus.name <- as.matrix(logs$gene)
locus.locus <- as.matrix(logs$locus)
logs <- t(logs[,-c(1,2,4)])

#GRowth
library(plyr)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
stage <- c("Egg", "Veliger", "Pediveliger", "Juvenile")
growth <- read.table("Growth.txt", sep="\t", header=T)
gm <- ddply(growth, .(DPF.1), numcolwise(mean))
gm.com <- rbind(c(0,5),				# size of the egg 5µm (?)
			apply(subset(gm[gm$DPF.1 %in% c(13,15),]), 2, FUN=mean),
			apply(subset(gm[gm$DPF.1 %in% c(18,20),]), 2, FUN=mean),
			gm[gm$DPF.1==25,]); gm.com$DPF.1 <- stage

#Performance
performance <- read.table("performance.txt", sep="\t", header=T)
pm <- ddply(performance, .(DPF.1), numcolwise(mean))
pm$DPF.1 <- stage
	
}
	
				# Calculation for One variable at a time
				#Linear regression of Gene1 on growth and performance
				{
				d <- logs[,2]; x <- gm.com[,2]; y <- pm[,2]
				summary(lm(d~x+y))
					# approach 1
					#calculate the standardized regression coefficients
					ds <- (d-mean(d))/sd(d)
					ys <- (y-mean(y))/sd(y)
					xs <- (x-mean(x))/sd(x)
					summary(lm(ds~xs+ys))
				# approach 2
				#Direct calculation of regression coefficients
				lm(d~x+y)$coefficients[2]*sd(x)/sd(d)
				lm(d~x+y)$coefficients[3]*sd(y)/sd(d)
				# d= 0.3949935y + 0.8571343x	R2= 0.675
}
## BATCH Calculation
#note d~x+y : coefficient[2]=x    & d~y+x : coefficient[2]=y
#The standardized regression coefficients for all variables (genes)
{
#Linear regression of ALL GENES (variables) versus growth and performance
x <- gm.com[,2]; y <- pm[,2]
#calculate ALL standardized regression coefficients (used for the regression biplot)
B <- NULL
for (j in 1:126) B <- rbind(B, lm(logs[,j]~y+x)$coefficients[2:3]*c(sd(y),sd(x))/sd(logs[,j]))

nrow(logs)
#anovaw2 j in 1:126
#anovaw-1 j in 1:106
#anovawNow j in 1:115
#anovawNob j in 1:103

str(summary(lm(logs[,2]~y+x)))
}

## Retrieve Pval, x, y
{
# Retrieve P-values (only the values)
P <- NULL
for (z in 1:126) P <- cbind(P,summary(lm(logs[,z]~y+x))$coefficients[2:3,4])
# Retrieve only significant P values
P <- as.data.frame(t(P))
Pv <- subset(P, x<='0.05' | y<='0.05')
Pv

# Retrieve significnat standardized coefficients (x, y)
rn <- rownames(Pv)				# significant lines in Pvalue
B.df <- as.data.frame(B); rnb <- rownames(B.df)				# rownames in B
lrn <- cbind(locus.name,rnb); lrn <- lrn[lrn[,2] %in% rn,]; lrn <- as.vector(lrn[,1])	#retrieve gene names
Lrn <- cbind(locus.locus,rnb); Lrn <- Lrn[Lrn[,2] %in% rn,]; Lrn <- as.vector(Lrn[,1])	#retrieve locus names
B.p <- as.data.frame(cbind(B, rnb))
B.s <- subset(B.p[B.p$rn %in% rn,1:2])
B.ss <- as.matrix((B.s))
B.ss <- as.matrix(cbind(as.numeric(B.ss[,1]),as.numeric(B.ss[,2])))
dimnames(B.ss) <- list(c(lrn),c("y", "x"))
B.ss
}

#Regression biplot
{
plot(xs, ys, xlab="x*(Growth)", ylab="y*(Performance)", type="n", asp=1, cex.axis=0.7,
xlim=c(-1.5,1.5), ylim=c(-1.5,1))
text(xs, ys, label=rownames(logs), col="grey", font=4)
points(B[,2:1], pch=3, cex=0.6,col="grey")
arrows(0,0,0.95*B.ss[,2],0.95*B.ss[,1],col="brown",angle=15,length=0.1)
text(B.ss[,2:1], label=rownames(Pv),adj=1.5, cex=0.7, col="red")
points(B.ss[,2:1], pch=21, type="p", bg="red")
}

# variance explained for each significant (rn) regression
{
var.fit <- NULL;var.expl <- NULL
for(j in unique(as.numeric(rn))) {var.fit <- cbind(var.fit,var(lm(logs[,j]~y+x)$fitted.values))
								print(var(lm(logs[,j]~y+x)$fitted.values)/var(logs[,j]))
								var.expl <- sum(var.expl + var(lm(logs[,j]~y+x)$fitted.values))
								}
var.expl/sum(apply(logs[,c(as.numeric(rn))],2,var))
}

## Retrieve the adjusted R^2 (to the p value) ~significnat only
{
R <- NULL
for (i in c(as.numeric(rn))) R <- rbind(R, summary(lm(logs[,i]~y+x))$adj.r.squared)
dimnames(R) <- list(c(rn),c("R2.adj"))
R
}

# Generate the final output
{
B.ss <- as.data.frame(cbind(Gene=rownames(B.ss),
							#locus=Lrn,					# If I want to add the locus ID
							y=round(as.numeric(B.ss[,1]),3),
							x=round(as.numeric(B.ss[,2]),3),
							variance=round(as.numeric(t(var.fit)),3),		#between y & x
							Rsq.adj=round(R,3)))
							}
B.ss
# Use stargaze (package) to draw LaTeX tables

#########################################################################
#								SECTION -6-								#
#########################################################################
=================
Chapter 3 : biplots GLM
=================
## both glms can't work with logarithmic numbers nor NAs
# Fatty acids dataset (no trocophore)

setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
acids <- read.table("ag.txt", sep="\t", header=T)
dataset <- acids[acids$Traitement %in% c("Cocktail"),]
library(plyr)
dataset <- ddply(dataset, .(Taille), numcolwise(mean))
dataset <- as.matrix(dataset[c(2,4,3,1),25:47])
nrow(dataset);ncol(dataset)


#GRowth
library(plyr)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
stage <- c("Egg", "Veliger", "Pediveliger", "Juvenile")
growth <- read.table("Growth.txt", sep="\t", header=T)
gm <- ddply(growth, .(DPF.1), numcolwise(mean))
gm.com <- rbind(c(0,5),				# size of the egg 5µm (?)
			apply(subset(gm[gm$DPF.1 %in% c(13,15),]), 2, FUN=mean),
			apply(subset(gm[gm$DPF.1 %in% c(18,20),]), 2, FUN=mean),
			gm[gm$DPF.1==25,]); gm.com$DPF.1 <- stage


#Performance
performance <- read.table("performance.txt", sep="\t", header=T)
pm <- ddply(performance, .(DPF.1), numcolwise(mean))
pm$DPF.1 <- stage
x <- gm.com[,2]; y <- pm[,2]
ys <- (y-mean(y))/sd(y)
xs <- (x-mean(x))/sd(x)


# Transformation (1) double Square root
# Double square root regressions
srt.dataset <- dataset[,1:23]^0.25
d.sqrt <- NULL;var.expl <-NULL
for(j in 1:23) {d.sqrt <- rbind(d.sqrt, lm(srt.dataset[,j]~ys+xs)$coefficients)
				var.expl <- sum(var.expl+var(lm(srt.dataset[,j]~ys+xs)$fitted.values))
				print(var(lm(srt.dataset[,j]~ys+xs)$fitted.values) / var(srt.dataset[,j]))
				}
var.expl / sum(apply(srt.dataset,2,var))
summary(lm(srt.dataset~ys + xs))				
				
plot(xs, ys, xlab = "x*(growth)", ylab = "y*(performance)", type = "n", asp = 1, cex.axis = 0.7)
text(xs, ys, labels = rownames(dataset), col = "forestgreen")
text(d.sqrt[,3:2], labels = colnames(dataset[,1:23]), col = "chocolate", font = 1)
for(j in 1:23) arrows(0, 0, 0.95*d.sqrt[j,3], 0.95*d.sqrt[j,2], col = "chocolate", angle = 15, length = 0.1)


# Trasformation (2) Poisson
#Performing poisson regression on all species

# Transformation (3) logit

#########################################################################
#								SECTION -7-								#
#########################################################################
=======================================
chapter 4: Multidimentional scaling MDS
=======================================
# Graphical representation of a set of objects based on their similarities
## Could be done on gene expressions. I should find how ?

#########################################################################
#								SECTION -8-								#
#########################################################################
===========================================================
chapter 5 : reduced dimensions (based on the rank of the matrix // SVD$d)
===========================================================
## SVD

## Work with 3 replicates Cocktail microarrays
## Genes
## Extract significant locus (anova) from original dataset meansW2
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
anova_SIG <- read.table("anovaW2.txt", sep="\t", header=F)
genes_SIG <- genes_coc[genes_coc$locus %in% anova_SIG$V2,]
logs <- as.matrix(t(genes_SIG[,4:18]))
dim(logs)	#15x126
}

# Singular value decomposition on genes
{
svd(logs)
#Generalized SVD on genes only
rowcoord <- svd(logs)$u %*% diag(sqrt(svd(logs)$d))
colcoord <- svd(logs)$v %*% diag(sqrt(svd(logs)$d))

plot(rbind(rowcoord,colcoord)[,1:15], type="n", xlab="", ylab="", asp=1, cex.axis=0.7)
abline(h=0,v=0,lty="dotted")
text(rowcoord[,1:15],labels=rownames(logs), col="forestgreen")
points(colcoord[,1:15], pch=3, cex=0.6)
}
#########################################################################
#								SECTION -9-								#
#########################################################################
======================
chapter 6 : PCA (unconstrained)
======================
## SVD & PCA


## NO pretransformation (no LRA)
## Work with 3 replicates COCktail
## Fatty acids dataset (no trocophore)
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
acids <- read.table("ag.txt", sep="\t", header=T)
fat <- acids[c(1:6,10:12,16:18),c(3,28:46)]
fat <- as.matrix(fat[,-1])
dim(fat)
}

## Centring & Weigthed SVD
## Matrix is multiplied by a constant
{
species <- NULL
for(j in 1:12) species <- rbind(species, (fat[j,]-(apply(fat,2,mean)))/(sqrt(nrow(fat)*ncol(fat))))
species
species.SVD <-svd(species)
}

## FORM BIPLOT
## Biplot (row=principal, col=std coordinates)
## Calculate coordinates
{
species.F <- sqrt(nrow(species))*species.SVD$u %*% diag(species.SVD$d)
species.G <- sqrt(ncol(species))*species.SVD$v
		## Contribution biplot - rescaling of the standard coordinates -
		#species.G <- ((ncol(species))^0.25)*species.SVD$v
plot(rbind(species.F[,1:2],species.G[,1:2]), type="n", asp=1,xlab="dim1",ylab="dim2", cex.axis=0.7)	#only 2 cols (even if there is more)
text(species.F[,1:12], labels=rep(c("E","V","P","J"),rep(3,4)),col="forestgreen", cex=0.8)
arrows(0,0,species.G[,1], species.G[,2],col="chocolate",length=0.1,angle=10)
text(1.2*species.G[,1], 1.2*species.G[,2],labels=colnames(fat),col="chocolate",font=4,cex=0.5)
}

## COVARIANCE BIPLOT (col=principal, row=std coordinates)
{
species.F <- sqrt(nrow(species))*species.SVD$u 
species.G <- sqrt(ncol(species))*species.SVD$v %*% diag(species.SVD$d)
plot(rbind(species.F[,1:2],species.G[,1:2]), type="n", asp=1,xlab="dim1",ylab="dim2", cex.axis=0.7)
text(species.F[,1:12], labels=rep(c("E","V","P","J"),rep(3,4)),col="forestgreen", cex=0.8)
arrows(0,0,species.G[,1], species.G[,2],col="chocolate",length=0.1,angle=10)
text(1.2*species.G[,1], 1.2*species.G[,2],labels=colnames(fat),col="chocolate",font=4,cex=0.5)
}

## eiganvalues
{
species.SVD$d^2
## sum of squares (SS) decomposed by the SVD = Total variance
## Total inertia
sum(species.SVD$d^2)*100
## SCREE PLOT
## Squared singular values are eiganvalues
eig.p <- 100*species.SVD$d / sum(species.SVD$d)
eig.p <- eig.p[seq(4,1)]
barplot(eig.p,horiz=T,cex.axis=0.7)
}
#########################################################################
#								SECTION -10-							#
#########################################################################
=======================================
chapter 7 : Logarithmic transformation	(unconstrained)
=======================================
## LRA & CA (SVD)


## Pretransofmration (LRA)
## Work with 3 replicates COCktail
## Fatty acids dataset (no trocophore, no sums)
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
acids <- read.table("ag.txt", sep="\t", header=T)
fat <- acids[c(1:6,10:12,16:18),45:47]			# keep only EFAs
}
		## Either ....
		fat <- acids[c(1:6,10:12,16:18),48:50]			# keep only total of the 3 classes of lipids


## Compute logs and decompose, next plot the results
## Summ of r and c
{
P <- fat/sum(fat)		
rm <- apply(P, 1, sum)			
cm <- apply(P, 2, sum)
## Logarithmic transformation
Y <- as.matrix(log(P))
	#na.omit(Y)	
## Weighted double centring
mc <- t(Y) %*% as.vector(rm)
Y <- Y - rep(1,nrow(P)) %*% t(mc)
mr <- Y %*% as.vector(cm)
Y <- Y - mr %*% t(rep(1,ncol(P)))
## calculating the centroid (length, width & height) = weighted average location of the matrix
Z <- diag(sqrt(rm)) %*% Y %*% diag(sqrt(cm))
svdZ <- svd(Z)				
## Calcuation of the principal and the standard coordinates
fat_F <- diag(1/sqrt(rm)) %*% svdZ$u[,1:2] %*% diag(svdZ$d[1:2])
fat_G <- diag(1/sqrt(cm)) %*% svdZ$v[,1:2]

plot(rbind(fat_F,fat_G/20), asp=1,type="n", xaxt="n", yaxt="n")
axis(1,col.axis="forestgreen", col.ticks="forestgreen")
axis(2,col.axis="forestgreen", col.ticks="forestgreen")
axis(3,col.axis="brown", col.ticks="brown")
axis(4,col.axis="brown", col.ticks="brown")
text(1.12*fat_F, labels=rep(c("E","V","P","J"),rep(3,4)), col="forestgreen")
points(fat_F, pch=3)
text(1.02*fat_G/20, labels=seq(ncol(fat)), col="brown")			# a little crowded
lines(fat_G/20, lty="dashed", col="chocolate", lwd=2)		# only for 3-5 variables
}
## Total variance
sum(svdZ$d^2)





