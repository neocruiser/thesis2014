#########################################################################
#								SECTION -1-								#
#########################################################################
## PCA section 1
## SVD & PCA

## GENES
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansW2.txt", sep="\t", header=T)
expressed <- as.matrix(genes[,4:18])
logs <- as.matrix(t(expressed))

## Centring & Weigthed SVD
## Matrix is multiplied by a constant
species <- NULL
for(j in 1:15) species <- rbind(species, (logs[j,]-(apply(logs,2,mean)))/(sqrt(nrow(logs)*ncol(logs))))
species.SVD <-svd(species)


## *********** COVARIANCE BIPLOT (col=principal, row=std coordinates) ****************
species.F <- sqrt(nrow(species))*species.SVD$u 
species.G <- sqrt(ncol(species))*species.SVD$v %*% diag(species.SVD$d)
plot(rbind(species.F[,1:2],species.G[,1:2]), type="n", asp=1,xlab="dim1",ylab="dim2", cex.axis=0.7)
abline(v=0, h=0, lty="dotted")
text(species.F[,1:15], labels=rep(c("E","T","V","P","J"),rep(3,5)),col="forestgreen", cex=0.6, font=4)
points(species.G[,1], species.G[,2], pch=20, col="grey", cex=0.6)

## Mark selected genes on the previous biplot
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
anova_SIG <- read.table("anovaW2.txt", sep="\t", header=F)
genes_SIG <- genes_coc[genes_coc$locus %in% anova_SIG$V2,]
rn <- genes_SIG$X
## Transcipts from ANOVA set
for(j in rn) points(species.G[j,1],species.G[j,2], pch=24, bg="orange", cex=0.5)

#points(species.G[14000,1],species.G[14000,2], pch=21, cex=0.7, bg="red")
#arrows(0,0,species.G[,1], species.G[,2],col="chocolate",length=0.1,angle=10)
#text(1.2*species.G[,1], 1.2*species.G[,2],labels=colnames(logs),col="chocolate",font=4,cex=0.5)
segments(species.G[,1:2])

## SCREE PLOT
## eiganvalues
species.SVD$d^2
## sum of squares (SS) decomposed by the SVD = Total variance accounted for by 2D solution
## Total inertia
sum(species.SVD$d^2)*100
## Average variance
sum(species.SVD$d^2)/(nrow(logs)-1)
## Total varaince explained by the two first axes (just sum up the first 2)
eig.p <- 100*species.SVD$d / sum(species.SVD$d)
eig.p <- eig.p[seq(15,1)]
barplot(eig.p,horiz=T,cex.axis=0.7)