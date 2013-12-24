#########################################################################
#								SECTION -1-								#
#########################################################################
## CCA between GENES (range) and FAT (constrained)	<- havent seen it anywhere else
## Greenakre

## Work with 3 replicates Cocktail microarrays
## Genes
## Extract significant locus (anova) from original dataset meansW2
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
anova_SIG <- read.table("anovaW2.txt", sep="\t", header=F)
genes_SIG <- genes_coc[genes_coc$locus %in% anova_SIG$V2,]
logs <- as.matrix(t(genes_SIG[,4:18]))
dim(logs)	#15x126
species <- decostand(logs, method="range")
species <- species[-c(4:6),]

setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
acids <- read.table("ag.txt", sep="\t", header=T)
fat <- as.matrix(acids[c(1:6,10:12,16:18),27:46])
dim(fat)
fat <- fat^0.25

library(vegan)
genes.cca <- cca(species,fat)

plot(genes.cca, display=c("lc","bp","sp"),type="n")
text(genes.cca, display="bp", col="chocolate", font=4, cex=0.6)
text(genes.cca, display="lc", col="blue", cex=0.6)

genes.cca.sp <- genes.cca$CCA$v.eig
genes.cca.spcon <- genes.cca$colsum * (genes.cca.sp[,1]^2 + genes.cca.sp[,2]^2) / sum(genes.cca$CCA$eig[1:2])
genes_names <- colnames(species)
genes_names[genes.cca.spcon<0.008] <- "."
	text(genes.cca, display="sp", labels=genes_names, col="forestgreen", cex=0.6)
		ordilabel(genes.cca, col = "black", dis = "sp", fill = "peachpuff", priority = genes_names, cex=0.4)
			orditorp(genes.cca, dis = "sp", prio = genes_names, pch = "+", pcol = "gray")

#########################################################################
#								SECTION -2-								#
#########################################################################
## Remove DOWnregulated genes
## keep only positive genes
## CCA genes + FAT
## strategy originated from watson 2006
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansW2.txt", sep="\t", header=T)
expressed <- as.matrix(subset(genes, o1>0 & o2>0 & o3>0 & v1>0 &v2>0 & v3>0 & p1>0 & p2 >0 & p3>0 &j1>0 &j2>0 &ja>0, select=(o1:ja)))
dim(expressed)
head(expressed)
species <- t(expressed[,-c(4:6)])
dim(species)

setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
acids <- read.table("ag.txt", sep="\t", header=T)
fat <- as.matrix(acids[c(1:6,10:12,16:18),27:46])
dim(fat)

library(vegan)
genes.cca <- cca(species,fat)

plot(genes.cca, display=c("lc","bp","sp"),type="n")
text(genes.cca, display="bp")
text(genes.cca, display="lc", col="red")

genes.cca.sp <- genes.cca$CCA$v.eig
genes.cca.spcon <- genes.cca$colsum * (genes.cca.sp[,1]^2 + genes.cca.sp[,2]^2) / sum(genes.cca$CCA$eig[1:2])
genes_names <- colnames(species)
genes_names[genes.cca.spcon<0.01] <- "+"
text(genes.cca, display="sp", labels=genes_names, col="forestgreen", cex=0.5)


#########################################################################
#								SECTION -3-								#
#########################################################################
## LRA with one env variable (CONSTRAINED CA)
## species = SVD, env= regression
## Greenacre

## Work with 3 replicates COCktail
## Fatty acids dataset (********no trocophore, no sums*********)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
acids <- read.table("ag.txt", sep="\t", header=T)
fat <- as.matrix(acids[c(1:6,10:12,16:18),27:46])

{
P <- fat/sum(fat)		
rm <- apply(P, 1, sum)			
cm <- apply(P, 2, sum)
Y <- as.matrix(log(P))
mc <- t(Y) %*% as.vector(rm)
Y <- Y - rep(1,nrow(P)) %*% t(mc)
mr <- Y %*% as.vector(cm)
Y <- Y - mr %*% t(rep(1,ncol(P)))
Z <- diag(sqrt(rm)) %*% Y %*% diag(sqrt(cm))
svdZ <- svd(Z)				
fat_F <- diag(1/sqrt(rm)) %*% svdZ$u
fat_G <- diag(1/sqrt(cm)) %*% svdZ$v %*% diag(svdZ$d)
}

## Weigths of larvae standardized on rm and regressed on the coordinates (cases)
## ************* Only continuous values ***********
weight <- acids[c(1:6,10:12,16:18),5]
weight.m <- cov.wt(as.matrix(weight), wt=rm)$center
weight.v <- cov.wt(as.matrix(weight),wt=rm)$cov
weight.s <- (weight-weight.m)/sqrt(weight.v)
lm(weight.s ~fat_F[,1]+fat_F[,2],weights=rm)$coefficients
summary(lm(weight.s ~fat_F[,1]+fat_F[,2],weights=rm))
## Constrained axis to be correlated w/ weight
W <- diag(sqrt(rm)) %*% weight.s %*% (1/(t(weight.s) %*% diag(rm) %*% weight.s)) %*% t(weight.s) %*% diag(sqrt(rm))
WZ <- W%*%Z
svdWZ <- svd(WZ)
## Orthogonal projection
WpZ <- Z - WZ
svdWpZ <- svd(WpZ)
## Principal (cols) & standard axes
fat_F[,1] <- diag(1/sqrt(rm)) %*% svdWZ$u[,1]
fat_G[,1] <- diag(1/sqrt(cm)) %*% svdWZ$v[,1] * svdWZ$d[1]
fat_F[,2] <- diag(1/sqrt(rm)) %*% svdWpZ$u[,1]
fat_G[,2] <- diag(1/sqrt(cm)) %*% svdWpZ$v[,1] * svdWpZ$d[1]
weight.c <- lm(weight.s~ fat_F[,1] + fat_F[,2], weights=rm)$coefficients


plot(rbind(fat_F,fat_G/2.5), asp=1,type="n", xaxt="n", yaxt="n")
axis(1,col.axis="forestgreen", col.ticks="forestgreen")
axis(2,col.axis="forestgreen", col.ticks="forestgreen")
axis(3,col.axis="brown", col.ticks="brown")
axis(4,col.axis="brown", col.ticks="brown")
text(1.12*fat_F, labels=rep(c("E","V","P","J"),rep(3,4)), col="forestgreen")
points(fat_F, pch=3)
text(fat_G/2.5, labels=colnames(fat), cex=0.6)

arrows(0,0,0.9*weight.c[1], 0.9*weight.c[2], lwd=1.5, length=0.1, angle=15)
text(weight.c[1], weight.c[2], "weight")

## Totl variance in the matrix
sum(Z*Z)
sum(WZ*WZ)
sum(WpZ*WpZ)
100*sum(WZ*WZ)/sum(Z*Z)
## The position of the observed percentage (over 10000 permutations)
set.seed(157)
custom_permut <- rep(0,10000)
total_var <- sum(Z*Z)
custom_permut <- NULL
## Start permutations
for(j in 1:10000) {
					weight.perm <- sample(weight.s)
					W <- diag(sqrt(rm)) %*% weight.perm %*% (1/(t(weight.perm) %*% diag(rm) %*% weight.perm)) %*% t(weight.perm) %*% diag(sqrt(rm))
					WZ <- W %*% Z
					custom_permut[j] <- 100*sum(WZ*WZ)/total_var
					}
pv <- which(sort(custom_permut)==custom_permut[1]); (10000-pv)/10000
