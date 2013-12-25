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
plotDensities(MA)
plotDensities(MA, col="black")

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
## Treatment main effect
cont.matrix <- makeContrasts(Veliger=(Vi-E)+(Vc-T),Pediveliger=(Pi-Vi)+(Pc-Vc),Juvenile=(Ji-Pi)+(Jc-Pc), levels=design)
## Stage of development main effect
cont.matrix <- makeContrasts(Veliger=(Vc+T)-(Vi-E),Pediveliger=(Pc+Vc)-(Pi-Vi),Juvenile=(Jc+Pc)-(Ji-Pi), levels=design)
## which genes respond differently over time between the two treatments ?
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

x11(title="PlotMA of MA linear model");plotMA(fit)
x <- topTable(fit2, adjust="BH",number=25,coef=1)
x
x<-topTable(fit2, adjust="BH",number=400)
x<-x[!duplicated(x$ProbeName),]
head(x);tail(x);dim(x)
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\datasets\\")
write.table(x$ProbeName,"ProbeNames_lmfit.txt",sep="\t",quote=F)

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

	# (20) edit output files to load them into MeV
	#remove quotation marks in notepad++ so as to get just 4 columns (variables) and 4 set of observations
	#remove numbering search: ^.*(CUST.*$) replace: \1

# ===========
## END
# ===========

#remove all objects
rm(list=ls())


# ========================================================================================
# Briefing before processing MeV output files
# ========================================================================================
# (step 1) = paste only CUSTs (MeV) into notepad++. Save File1
# (step 2) = type I file = retrieve Locus ID of CUST. Save File2
# (step 3) = type II file = retrieve nucleotide sequence for significant MeV spots. Save File3
# (step 4) = type III file = remove sequence, keep locus|GO ID|name. Save File4
# (additional step 5) = paste ANOVA data (MeV). Save File0
# (additional step 6) = import File0 and File2 in excel


## START
# ========================================================================================
# Reason: create a type I file.
# ========================================================================================
# Retrieve "locus names" to "CUSTs" selected with MeV
# (1) add "cust" column in cust.txt files

library(limma)
setwd("C:/Dropbox/Workshop2013/Work/Microarrays/10- MeV/Gene Selection/")
id.t <- readTargets("id.txt")         #from bacterial spot removal
# (2) change only file names then repeat the following codes
CUST <- readTargets("anova1.txt")   #input file (name); CUSTs generated w/ MeV
CUST.t <- CUST[!duplicated(CUST$cust),]    #remove duplicates, hence creating a data frame of "characters" which could be read by %in%

locus <- id.t[id.t$cust %in% c(CUST.t),]

x <- as.data.frame(locus)             #convert to "data frame" because "%in% locus" is a "character matrix"
locus.t <- x[!duplicated(x$locus),]         #only data frames could be treated for duplication
write.table(locus.t, "anova2.txt", sep="\t")
rm(CUST, CUST.t, locus, x, locus.t)
## END


## START
# ========================================================================================
# Reason: Create type II files. w/o SAM. w/o CUST. Retrieve "Nucleotide sequence"
# ========================================================================================
# Find nucleotide sequences from differentially expressed genes (CUSTs from MeV output)
# edit the "export_seq_blast15457.fasta" (build a data frame that can be imported in R, reason to find sequences specific to MeV CUSTs)
# (0) check that there is only 15457 occurences being replaced
# (1) notepad++ search: (^>.*$)\r\n([^>]*$) and replace: \1\t\2 (in regular expression mode)
# (2) notepad++ search: \| and replace: \t\|
# (3) notepad++ search: ^> and replace: >\t
# (4) add first line (x locus description sequence)
# (5) save as seperate file "exportID.txt"
# (6) when finished, reverse repeat steps (3), (2), (1)
# (7) notepad++ search: " and replace: 'blank'
# (8) notepad++ search: ^>\t and replace: ^>
# (9) notepad++ search: \t\| and replace: \|
# (10) notepad++ search: (^>.*)\t and replace: \1\r\n
# (hint) Procedure done once per dataset. Meaning exportID is for all sequences

library(limma)
setwd("C:/Dropbox/Workshop2013/Work/Microarrays/10- MeV/Gene Selection/")
ID <- readTargets("exportID.txt")        # containing the 15k annotated Blast2GO sequences
# (7) change only file names then repeat the following codes
exp <- read.table("anova2.txt", header=TRUE, sep="")      # containing the MeV exported experimental locuses (SAM statistics)
locus <- as.character(exp[,2])
seq <- ID[ID$locus %in% c(locus),]
write.table(seq, "anova2.txt", sep="\t")
rm(exp, seq, locus)
## END

## START
# ========================================================================================
# Reason: create type III files
# ========================================================================================
# Traitement des séquences après les avoirs indexées puis extraits, create a tab delimited files
# (O) copy AnnotSeq.txt files into another directory
# (update 1) notepad++ search: .*?> and replace: >
# (update 2) notepad++ search: " and replace: blank
# (update 3) notepad++ search: >\t and replace: >
# (update 4) notepad++ search: \t| and replace: |
# (update 5) notepad++ search: \t and replace: \r\n
# (1) notepad++ search: ^[^>]*$ and replace: 'blank' (removes all nucleotide sequences)
# (2) notepad++ search: \r\n> and replace: 'blank' (removes all blank lines)
# (4) notepad++ search: | and replace: \t (adds a tab between variables)
# (not necessary anymore) notepad++ search: (?<=\t)(.*$) and replace: "\1"
# (6) notepad++ add header names (oeuf  trocophor)(oeuf veliger)(oeuf pediveliger) etc.
# (hint) All files must be sep=""
# (tip) R will take every word seperated by a space as a new variable then prompt an error, this is why you should add quotation mark to every observation
# (7) make sure you have:
# (A) File set for CUST only
# (B) File set for CUST vs locus
# (C) File set for locus vs descriptin (w/o sequence)
# (D) File set for SAM MeV data
# (E) File set for CUST - locus - description - sequence - SAM
# (hint) Merge B - C - D into one excel file
# (tip) remove space between 2 words, if so, add quotation marks (eg, "hello world")
# (tip) type I files contains 3 variables, 1 of which should be deleted.
# (1) notepad++ search: ^.*\t(CUST.*?$) and replace: \1
# type I files == CompleteCust.txt
# type II files == CompleteLocus.txt
# type III files == CompleteSAM.txt

library(limma)
setwd("C:/Dropbox/Workshop2013/Work/Microarrays/10- MeV/Gene Selection/")
data.path <- "C:/Dropbox/Workshop2013/Work/Microarrays/10- MeV/Gene Selection/"
datafiles <- list.files(path=".", pattern="locus_.*?.txt$")
data <- do.call(rbind, lapply(datafiles, read.table))
write.table(data, "CompleteCust.txt", sep="\t")
rm(datafiles, data)
