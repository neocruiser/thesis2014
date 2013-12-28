=========================================================================
#		SECTION --								#
=========================================================================
# ===================================================================================================
# Preprocessing of time-course design of Mussel Larvae Mytilus edulis. Ph.D 2014. Agilent microarray
# ===================================================================================================

setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\datasets\\")
load("paper3.Rdata", .GlobalEnv)
lsos(pat="")

library(limma)
library(statmod)
# (3) CREATE LIMMA TARGETS
set.seed(2013)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
targets <- readTargets("targets4.txt")
RG <- read.maimages(targets, source="agilent")
# (5) BACKGROUND CORRECTION
RG_min <- backgroundCorrect(RG, method="minimum")	## PAPER 2
RG_normexp60 <- backgroundCorrect(RG, method="normexp", offset=60, normexp.method="mle")	# PAPER 3e

n=sample(1:nrow(targets),6)		## number of arrays in target file
## PLOT (1)
lsos(pat="")
v0 = RG			## set v0 to any RG matrix
windows(width=10, height=15, title="RG normexp10", restoreConsole=TRUE)
par(mfrow=c(3,2))
for(i in n){
	plotMA(v0, array=i)
abline(h=0, col="red", lty=2, lwd=1)
}

## PLOT (2)
x11(); plotDensities(v0)
graphics.off()


# =========================================================================================================
# (6) setting spot quality weights (flagged out spots will not influence the normalization of other spots)
# =========================================================================================================
# (7) normalize w/n arrays then remove probes
# (tip) reference to "regular expression" on wikipedia
# (A) notepad++ search: (TERA).*?(?=\t) (regular exp mode) and replace: titration (action to replace TERA0123 string with titration)
# (hint) '.*?' means any character after TERA. all of them including '\t'
# (hint) '(?=\t)' means stop at the end and don't include '\t' tabulation
# (B) (?<=\t).*(Corner).*?(?=\t) for DarkCorner search/replace w/ Corner or directly by DarkCorner & GE_BrightCorner
# (C) repeat step (A) while replacing TERA w/ E1A, (+)E1A, ERCC, ETG, DCP
# (tip) setup on all arrays in the target.txt
# (hint) upweight non differentially expressed spots between arrays (p29 limma userguide)

w <- modifyWeights(array(1,dim(RG)),
                   RG$genes$SystematicName,
				   c("3xSLv1", "E1A", "ERCC", "ETG","titration"),
				   c(-1,-1,0,0,2))

lsos(pat="")
v1 = RG_normexp60
MA_normexp_NormW <- normalizeWithinArrays(v1,method='loess')	# PAPER 3

MA_normexp_NormWrobsp <- normalizeWithinArrays(v1, method="robustspline")		## just for experimentation
## method = control w weighting automatically
ctrl <- RG$genes$SystematicName=="titration"; summary(ctrl)
MA_normexp_NormWctrl <- normalizeWithinArrays(v1, method="control",controlspots=ctrl)		## just for experimentation
MA_normexp_NormWw <- normalizeWithinArrays(v1, weights=w)	## using weightse

## PLOTTING (1)
lsos(pat="")
v0 = MA_normexp_NormW			## set v0 to any MA matrix
windows(width=10, height=15, title="RG plain|nomW()")
par(mfrow=c(3,2))
for(i in n){
	plotMA(v0, array=i)
	abline(h=0, col="red", lty=2, lwd=1)
}


				# (10) array quality control
				# (hint) setting empirical array weights (p64 limma userguide)
				arrayw <- arrayWeights(MA)
				arrayw
				barplot(arrayw, width=c(0.2,0.2), xlim=c(0,7), xlab="Array", ylab="Weight", col="white", las=2)
				abline(h=1, col="black", lty=2, lwd=1)
				# (11) ebayes after setting quality array weights
				fitw <- lmFit(MA, weights=arrayw)
				fitw <- eBayes(fitw)
				topTable(fitw, number=30, adjust="BH", genelist=fitw$genes$SystematicName, p.value=0.05)
				volcanoplot(fitw, highlight=15)
				# (hint) le boxplot prendra en compte le status 0, +1 et -1 des probes
				fit <- lmFit(MA)
				fit <- eBayes(fit)
				boxplot(fit$t~MA$genes$ControlType,
						at=3:5-0.2, col=5, boxwex=0.4, xlab="control type",
						ylab="moderated t-statistics", pch=".", ylim=c(-20, 20), medlwd=1)
				boxplot(fitw$t~MA$genes$ControlType,
						at=3:5+0.2, col=6, boxwex=0.4, add=TRUE, yaxt="n", medlwd=1, pch=".")
				abline(h=0, col="black", lty=2, lwd=1)
				legend("topleft", legend=c("Equal weights", "Array weights"), fill=c(5,6), cex=0.8)
				vennDiagram(fit, include="up")


# (12) normalize between arrays
# not appropriate here see page 83 of limma userguide
# (hint) scale normalization to equalize the median absolute value for each array
lsos(pat="")
v1 = MA_normexp_NormW
MA <- normalizeBetweenArrays(v1, method="quantile")	## PAPER 3 (paper 2 was wo normalization Between)

MA <- normalizeBetweenArrays(v1, method="Gquantile")
MA <- normalizeBetweenArrays(v1, method="Tquantile", targets=targets)

## PLOT (1)
x11();plotPrintTipLoess(MA)
x11();plotDensities(MA)
#plotDensities(MA, col="black")

## PLOT (2)
x11()
cols <- colnames(MA$M)
cols[1:3] <- "red"
cols[4:6] <- "blue"
cols[7:9] <- "violet"
cols[10:12] <- "yellow"
cols[13:15] <-  "grey"
cols[16:18] <- "violet"
cols[19:21] <- "yellow"
cols[22:23] <- "grey"
boxplot(MA$M~col(MA$M),names=colnames(MA$M),col=cols,xlab="Mussel Larvae",ylab="M-values")

## PLOT (3)
windows(width=10, height=15, title="RG plain|nomW()")
par(mfrow=c(3,2))
lsos(pat="")
v0 = MA			## set v0 to any MA matrix
n=nrow(targets)
for(i in 1:n){
plotMA(v0, array=i, las=1)
abline(h=0, col="chocolate", lty=2, lwd=1)
x <- MA[MA[[2]]$SystematicName=="titration",] ; points(x[[6]], x[[5]], col="brown", pch=".")
x <- MA[MA[[2]]$SystematicName=="ERCC",] ; points(x[[6]], x[[5]], col="blue", pch=".")
x <- MA[MA[[2]]$SystematicName=="3xSLv1",] ; points(x[[6]], x[[5]], col="yellow", pch=".")
x <- MA[MA[[2]]$SystematicName=="ETG",] ; points(x[[6]], x[[5]], col="pink", pch=".")
#x <- MA[MA[[2]]$SystematicName=="E1A",] ; points(x[[6]], x[[5]], col="grey", pch=".")
x <- MA[MA[[2]]$SystematicName=="corner",] ; points(x[[6]], x[[5]], col="green", pch=".")
x <- MA[MA[[2]]$SystematicName=="DCP",] ; points(x[[6]], x[[5]], col="orange", pch=".")
}
dim(MA)

# (14) remove control probes spike-ins (+1 and-1) // E1A,Corner,ERCC,ETG,DCP,3xSLv1 (1319 spot)
dim(MA)
MA <- MA[MA$genes$ControlType=="0", ]
dim(MA)
# (15) remove all non-annotated sequences (27744 spot)
MA <- MA[!MA$genes$SystematicName=="---NA---", ]
dim(MA)
# (16) remove TERA sequences (1400 spot)
MA <- MA[!MA$genes$SystematicName=="titration", ]
dim(MA)

## PLOT
windows(width=10, height=15, title="RG plain|nomW()")
par(mfrow=c(3,2))
for(i in n){
	plotMA(MA, array=i)
	abline(h=0, col="red", lty=2, lwd=1)
}
dim(MA)

## SAVE
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
lsos(pat="")
save(list=ls(pattern=""),file="paper3.Rdata")

=========================================================================
#		SECTION --						#
=========================================================================

## experimental with Lmfit function (If not interested continue at section 5)
## (17) design matrix + linear model (page 53 of limma userguide)
## first design (must run contrast.fit)
lev <- c("E","T","Vc","Pc","Jc","Vi","Pi","Ji")
f <- factor(targets$Cy5)
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
design

## choose one contrast matrix ...
## Treatment main effect (Contrast 1)
cont.matrix <- makeContrasts(Veliger=(Vi-E)+(Vc-T),Pediveliger=(Pi-Vi)+(Pc-Vc),Juvenile=(Ji-Pi)+(Jc-Pc), levels=design)
## Stage of development main effect (Contrast 2)
cont.matrix <- makeContrasts(Veliger=(Vc+T)-(Vi-E),Pediveliger=(Pc+Vc)-(Pi-Vi),Juvenile=(Jc+Pc)-(Ji-Pi), levels=design)
## which genes respond differently over time between the two treatments ? (Contrast 3)
cont.matrix <- makeContrasts(Veliger=(Vi-E)-(Vc-T),Pediveliger=(Pi-Vi)-(Pc-Vc),Juvenile=(Ji-Pi)-(Jc-Pc), levels=design)
## which genes respond at either stage of the COCKTAIL treatment ?
## which genes respond at either stage of the ISO treatment ?
## which genes differently expressed relatively to EGG phase


## Group mean parametrization
fit <- lmFit(MA,design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
## PLOT VennDiagram
coefTests <- decideTests(fit2)
vennDiagram(coefTests)

x11(title="PlotMA of MA linear model");plotMA(fit)
x <- topTable(fit2, adjust="BH",number=25,coef=1)	## coef is dependent on the ncol(contrast.matrix)
x
x<-topTable(fit2, adjust="BH",number=400)
x<-x[!duplicated(x$ProbeName),]
head(x);tail(x);dim(x)

## EXTRACT 2FOLD GENES
## (hint) results is the preprocessed dataframe in dataset.R
names(results)
logsLM <- results[results$cust %in% x$ProbeName,]
logsLM <- as.matrix(logsLM)
logsLM <- logsLM[! duplicated(logsLM[,(dim(logsLM)[2]-1)]),]
cat("\n","All Toptable genes:",dim(logsLM)[1],"\n")
logs_2fold <- NULL
for(i in 2:dim(logsLM)[2]){
logs_2fold <- rbind(logs_2fold, logsLM[logsLM[,i] >'0.9',])
}
logs_2fold <- logs_2fold[! duplicated(logs_2fold[,(dim(logs_2fold)[2]-1)]),]
cat("\n","Only 2fold minimum expression genes:",dim(logs_2fold)[1],"\n")
sm <- sample(1:dim(logs_2fold)[1],5)
logs_2fold[sm,]
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\datasets\\")
write.table(logsLM, "Contrast3.txt", sep="\t", quote=F)

## VENN FOR CONTRAST TESTING
require(gplots)
contrast1 <- as.vector(logsLM[,24])
contrast2 <- as.vector(logsLM[,24])
contrast3 <- as.vector(logsLM[,24])
input <- list(Contrast1=contrast1,Contrast2=contrast2,Contrast3=contrast3)
venn(input)

# (19) add missing columns and build table of all selected genes
results1 <- MA$genes$ProbeName
results2 <- MA$M
results3 <- cbind(results1, results2)       #merge columns
colnames(results3) <- MA$genes$names
results <- data.frame(results3)
names(results) <- c("cust", "E1", "E2", "E3","T1","T2","T3","Vc1","Vc2","Vc3","Pc1","Pc2","Pc3","Jc1","Jc2","Jc3","Vi1","Vi2","Vi3","Pi1","Pi2","Pi3","Ji1")	## PAPER 3
#names(results) <- c("cust", "E1", "E2", "E3","T1","T2","T3","Vc1","Vc2","Vc3","Pc1","Pc2","Pc3","Jc1","Jc2","Jc3")	## PAPER 2
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\datasets\\")
write.table(results, "resultsCoc_Iso.txt", sep="\t", quote=F)

## PLOT histograms
## (hint) results is the preprocessed dataframe in dataset.R
require(lattice)
par(mfrow=c(1,3))
histogram(~E1+E2+E3+T1+T2+T3, data=results)
histogram(~Vc1+Vc2+Vc3+Vi1+Vi2+Vi3, data=results)
histogram(~Pc1+Pc2+Pc3+Pi1+Pi2+Pi3, data=results)
histogram(~Jc1+Jc2+Jc3+Ji1, data=results)

## PLOT (2) -- experimental
plot(Vc1~E1, data=results, pch=19, cex=0.5)
points(results$Vc1, results$T1, pch=19, cex=0.5, col="red")


# ===========
## END
# ===========

#remove all objects
rm(list=ls())
