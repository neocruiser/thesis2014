# Paul Pavlidis (c) 2003-2010

#########################################################################
#								SECTION -1-								#
#########################################################################

## one way anova
## READ data
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T, row.names=3)
sdata <- as.matrix(genes[,3:17])

	#species <- gl(5,3,15, label=c("E", "T", "V", "P", "J"))
	#sampl <- gl(5,3,15, label=paste("s", as.character=seq(1,5),sep=""))
	## species must be balanced : same number of terms in each groups
species <- gl(2,9,15, label=c("LA","PL"))	## I can replace label by levels()
sampl <- gl(5,3,15, label=c("E", "T", "V", "P", "J"))

# define ANOVA function
aof <- function(x) { 
  m<-data.frame(species, sampl, x); 
  anova(aov(x ~ sampl + species  + sampl*species, m))}
anovaresults <- apply(sdata, 1, aof)
pvalues<-data.frame( lapply(anovaresults, function(x) { x["Pr(>F)"][1:3,] }) )
str(head(pvalues))
setwd("C:/Dropbox/Workshop2013/Work/R/ANOVA/")
write.table(t(pvalues), file="anova-results3.txt", quote=F, sep="\t")
pvalues[,1:5]
}

# Benjamini_Hochberg FDR analysis
{
bh <- function(x, fdr) {
  thresh <- F;
  crit<-0;
  len<-length(x)
  answer <- array(len);
  first <- T;
  for(i in c(len:0)) {
    crit<-fdr*i/len;
    if (x[i] < crit || thresh == TRUE) {
      answer[i]<-T
      thresh <- T
      if (first) {
        cat(i ,"genes selected at FDR =", fdr ,"\n")
        first = F;
      }
    } else {
      answer[i]<-F
    }
  }
  answer
}
regionp <- sort(t(pvalues[1,]))
}
## *********HOLD********** set up a threshold to view the FDR in this anova.
## only changing the anova will change the FDR
fdr.result <- bh(regionp, 0.05681)
## write (optional)
{
bhthresh <- cbind(regionp, fdr.result)
write.table(bhthresh, "bhthresh2.txt", sep='\t', quote=F)
}

## select high pvalues 
## ****HOLD*********SET pVALUES (15257*0.0007/nrow(reg.hi.p))
reg.hi.p <-t(data.frame(pvalues[1, pvalues[1,] < 0.0007]))
## Run
{
cat(nrow(reg.hi.p) ," ... number of significant genes","\n")
reg.hi.pdata <- sdata[ row.names(reg.hi.p), ]
}

## ***********HOLD
## set template match
cbtempl<-c(0,0,0,0,0,0,1,1,1,0,0,0,0,0,0)
{
## Set gene by hand
cor.test((reg.hi.pdata[5,]),cbtempl)

template.match <- function(x, template) {
  k<-cor.test(x,template)
  k$p.value
}
cbtempl.results <- apply(reg.hi.pdata, 1, template.match, cbtempl)
temp.sig <- as.matrix(sort(cbtempl.results))

## Extract all logs for the template matches
NN <- as.matrix(temp.sig[temp.sig[,1]<0.05,1])
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansLmFitBgmin.txt", sep="\t", header=T, row.names=3)
logs_temp <- genes[rownames(genes) %in% rownames(NN),]
logs_temp <- t(logs_temp[,-c(1:2)])
cat("\n","\n",nrow(NN),"... number of significant genes (@ .05) that match the setup template","\n")
}
## write (optional)
setwd("C:/Dropbox/Workshop2013/Work/R/ANOVA/")
write.table(temp.sig, "template-results.txt", quote=F, sep="\t")
write.table(logs_temp, "meansLmFitBgmin_VEL.txt", quote=F, sep="\t")

## ***********HOLD
## set sampl= (hint:same as template)
t.test(t(reg.hi.pdata[1,sampl!="E"]), t(reg.hi.pdata[1,sampl=="E"]), var.equal=T)
## write (optional)
{
setwd("C:/Dropbox/Workshop2013/Work/R/anova2/")
write.table(reg.hi.pdata[order(cbtempl.results),], "cbtempl.pdata.txt", quote=F, sep='\t')
}



## *************STOP***************
## PREPARE CLUSTERING
## Averaging Samples (by col)
## Extract logs
x <- data.frame(Y=paste(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))),t(reg.hi.pdata))	## pvalues ALL
## either one
x <- data.frame(Y=paste(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))),logs_temp)		## pvalues TEMPLATE
{
library(plyr)
y <- ddply(x, .(Y), numcolwise(mean))
## Final dataset (average col)
pdata_hi <- t(y[,-1])
colnames(pdata_hi) <- c("Egg", "Trocophore", "Veliger", "Pediveliger", "Juvenile")
head(pdata_hi)
}

## create a minimal ExpressionSet & run MFUZZ
{
library(Biobase)
minimalSet_ns <- ExpressionSet(assayData=pdata_hi)
## Soft clustering and visualisation
library(Mfuzz)
## Standardization (relative expression vs absolute ampliture)
minimalSet <- standardise(minimalSet_ns)
## Fuzzifier m testing
mestimate(minimalSet)
}
## RUN soft clustering
cl <- mfuzz(minimalSet,c=6,m=1.5)
mfuzz.plot(minimalSet,cl=cl,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))


## coupling between clusters (high values between clusters = smilar pattern)
{
O <- overlap(cl)
over_plot <- overlap.plot(cl, over=O, thres=0.05)
cat("\n","\n","If coupling is large, cluster patterns are similar","\n")
over_plot
## Analysis of global clustering structures
cl2 <- mfuzz(minimalSet,c=6,m=1.5)
mfuzz.plot(minimalSet,cl=cl2,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
O2 <- overlap(cl2)
over_plot <- overlap.plot(cl2, over=O2,P=over_plot, thres=0.05)
}
## ADDITIONAL OPTIONS (EXTRACT CLUSTERS)
## See available components & extract GENES
{
cl
## Plot the center of the third cluser
plot(cl[[1]][3,],type="l",ylab="Expression")
## Extract the first 10 membership of the first cluster
cl[[4]][1:10,1]
## Extract genes forming the alpha cores of soft clusters
acore_list <- acore(minimalSet, cl=cl, min.acore=0.7)
str(acore_list)
cluster1 <- as.matrix(acore_list[[1]]$NAME)
cluster2 <- as.matrix(acore_list[[2]]$NAME)
cluster3 <- as.matrix(acore_list[[3]]$NAME)

}

#########################################################################
#								SECTION -2-								#
#########################################################################
# Commands demonstrating the use of R for statistical analysis of gene
# expression microarray data. This may be freely distributed, modified
# and used, with attribution.

# This code has been tested for R 1.6.2 under Windows 2000 and Unix
# (Solaris, Linux), and updated for R 2.11 (November 2010)

# Note that to avoid problems with invalid R variable names, in the
# input data file all "_" have been replaced with ".". Thus in this
# demo probe sets are named XXXXX.at instead of the usual Affymetrix
# (tm) format XXXXX_at.

# Under Windows, use the "File...change dir" menu command to go to the
# directory where the sample file can be found.
# Read in the data, set up variables
{
sdata<-read.table("sandberg-sampledata.txt", header=T, row.names=1)
strain <- gl(2,12,24, label=c("129","bl6"))
region <- gl(6,2,24, label=c("ag", "cb", "cx", "ec", "hp", "mb"))

# define ANOVA function
aof <- function(x) { 
  m<-data.frame(strain,region, x); 
  anova(aov(x ~ strain + region + strain*region, m))
}

# apply analysis to the data and get the pvalues.
anovaresults <- apply(sdata, 1, aof)

# Note the the following line works on the test data set, but fails
# with a stack overflow on large data sets, due to an apparent bug in
# R.
pvalues<-data.frame( lapply(anovaresults, function(x) { x["Pr(>F)"][1:3,] }) )

# (Begin bug work-around section)
# Workaround 1: The following can be used to replace the previous line.
# Suggested by Tom Blackwell, U Michigan.
{
temp.1 <- unlist(lapply(anovaresults, function(x) { x["Pr(>F)"][1:3,] }))
temp.2 <- matrix(temp.1, length(anovaresults), 3, byrow=T)
dimnames(temp.2) <- list(names(anovaresults), dimnames(anovaresults[[1]])[[1]][1:3])
pvalues <- data.frame(t(temp.2))

# Workaround 2 (also suggested by T. Blackwell). Replace the aof
# function with the following:
aof <- function(x) { 
  m<-data.frame(strain,region, x); 
  temp<-anova(aov(x ~ strain + region + strain*region, m))
  temp["Pr(>F)"][1:3,]
}
# And then do the following:
pvalues<-data.frame(apply(sdata, 1, aof));
# Then proceed as described below.  The drawback of the preceding is
# that other data in the aov objects (e.g., F statistics) are thrown
# away.
}
# (End bug work-around section)

# Note that the results of all calls to write.table do not yield
# correct headers for use in most other programs. This is because no
# label is placed in the "upper left hand corner" of the matrix. When
# programs like Excel read in such files, the line containing the
# column labels is shifted to the left. This can be repaired by using
# Excel or even a word processor to insert a heading such as "Probe"
# at the start of the first row.
setwd("C:/Dropbox/Workshop2013/Work/R/anova2/")
write.table(t(pvalues), file="anova-results.txt", quote=F, sep='\t')


# Get the genes with good region effect pvalues. # Note: removed sorting 
reg.hi.p <-t(data.frame(pvalues[2, pvalues[1,] < 0.001 & pvalues[2,] > 0.1]))
reg.hi.pdata <- sdata[ row.names(reg.hi.p), ]

# Template matching
cbtempl<-c(0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0)
cor.test(t(reg.hi.pdata[14,]),cbtempl)

# Run template on all the high region effect genes
template.match <- function(x, template) {
  k<-cor.test(x,template)
  k$p.value
}
cbtempl.results <- apply(reg.hi.pdata, 1, template.match, cbtempl)
write.table(cbtempl.results, file="cbtempl-results.txt", quote=F, sep='\t')

# t-test demonstration.
t.test(t(reg.hi.pdata[1,region!="cb"]), t(reg.hi.pdata[1,region=="cb"]), var.equal=T)

# output the region genes in order of cerebellum match.
write.table(reg.hi.pdata[order(cbtempl.results),], "cbtempl.pdata.txt", quote=F, sep='\t')

# Benjamini_Hochberg FDR analysis
bh <- function(x, fdr) {
  thresh <- F;
  crit<-0;
  len<-length(x)
  answer <- array(len);
  first <- T;
  for(i in c(len:0)) {
    crit<-fdr*i/len;
    if (x[i] < crit || thresh == T) {
      answer[i]<-T
      thresh <- T
      if (first) {
        cat(i ,"genes selected at FDR =", fdr ,"\n")
        first = F;
      }
    } else {
      answer[i]<-F
    }
  }
  answer
}

regionp <- sort(t(pvalues[2,])) # the p values must be sorted in increasing order!
fdr.result <- bh(regionp, 0.05) # reports that 192 genes are selected
bhthresh <- cbind(regionp, fdr.result)
write.table(bhthresh, "bhthresh.txt", sep='\t', quote=F) # print to a file.

# Figure 3a: Make a graph and save as postscript (these commands are not shown in the paper)
barplot(as.matrix(reg.hi.pdata[order(cbtempl.results),][1,]), names.arg=as.vector(sampl), xlab="Sample", ylab="Expression", col="grey", space=0)
dev.copy2eps(file="example.ps")

# Figure 3b, users of matrix2png can run the shell command (unix or
# macosX) or use the matrix2png web interface at
# http://microarray.cpmc.columbia.edu/matrix2png/ on
# cbtempl.pdata.txt, choosing the 'normalize' options, choosing color
# map 1 and setting the range to be -2 to 2. Note that to be used in
# matrix2png, the header of cbtempl.pdata.txt must be corrected as
# described above, and saved as text. For publication, the resulting
# figure was relabeled in Adobe Illustrator to make figure 3b. In R
# for unix, assuming matrix2png is installed on your system:
system("matrix2png -data cbtempl.pdata.txt -map 1 -zrcs -size 8:18 -numr 10 -range -2:2 > cbtempl.pdata.png")
}


#########################################################################
#								SECTION -3-								#
#########################################################################
## Use a pattern from subset C SAM
## pattern discovery with quantitative trait
## pattern (0,0,0,0,1)

## Extract selected locuses
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
## UP genes
siggenes_up <- read.table("quanti_siggenes_up.txt", sep="\t", header=T)
locus_up <- siggenes_up[,c(2)]
			## LOW genes
			siggenes_lo <- read.table("quanti_siggenes_lo.txt", sep="\t", header=T)
			locus_lo <- siggenes_lo[,c(3,2)]

			sam_pat <- rbind(locus_up[,2],locus_lo[,2])

setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansW2.txt", sep="\t", header=T, row.names=3)
sdata <- as.matrix(genes[,3:17])
## dataset for anova computation
sdata <- sdata[rownames(sdata) %in% locus_up,]

