#########################################################################
#								SECTION -1-								#
#########################################################################*
## DENDROGRAMS FOR EVERY CLUSTER (OBTAINED WITH MFUZZ)
## PAPER 2
## set color
?RColorBrewer
library(RColorBrewer)
display.brewer.pal(9, "YlGnBu")
palett<-brewer.pal(9,"Greens")

## Load SoftClustering Mfuzz package --- 12 clusters
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
load("clusters12.Rdata", .GlobalEnv)
lsos(pat="cluster.*")


k=cluster3
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes <- read.table("meansW2.txt", sep="\t", header=T)
genes.c <- genes[genes$locus %in% k[,1],]
names(genes.c)
rownames(genes.c) <- paste("T",rep(1:length(genes.c[,2])),genes.c[,2])					## change names
genes.c <- genes.c[,-c(1:3)]
require(plyr)
samples <- as.numeric(gl(5,3,15, label=seq(1:5)))
genes.cd  <- data.frame(samples=samples, t(genes.c))
genes.cd <- ddply(genes.cd, .(samples), numcolwise(mean))
genes.cd <- t(genes.cd)

genes.co <- as.matrix(genes.cd[order(genes.cd[,1]),])
genes.co <- genes.co[! rownames(genes.co)=="samples",]
windows(width=30,height=30)
heatmap.g <- heatmap(genes.co, Rowv=NA, Colv=NA, col=palett, margins=c(1,40), scale="row", labCol=c("E","T","V","P","J"))

#########################################################################
#					SECTION -2-								#
#########################################################################
## HIERARCHICAL AND BOOTSTRAP ANALYSIS
## (help) source : http://tinyurl.com/lqs76tf
require(pvclust)
require(gplots)

## load sets
setwd("C:/Dropbox/Workshop2013/Work/R/ANOVA/")
load("anova.Rdata", .GlobalEnv)
lsos(pat="nova")
## set NxN
NxN = setAnova

myanova <- t(NxN)
mydataanova <- t(scale(t(myanova)))
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
hra <- hclust(as.dist(1-cor(t(mydataanova), method="pearson")), method="complete") ## ROWS (genes)
hca <- hclust(as.dist(1-cor(mydataanova, method="pearson")), method="complete")	## COL (samples)
heatmap(myanova, Rowv=as.dendrogram(hra), Colv=as.dendrogram(hca), col=my.colorFct(), scale="row")

## CUT THE TREE
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
mycl <- cutree(hra, h=max(hra$height)/2)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(myanova, Rowv=as.dendrogram(hra), Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", RowSideColors=mycolhc)

## BOOTSTRAPING
pvAnova <- pvclust(scale(t(myanova)), method.dist="correlation", method.hclust="complete", nboot=2000)
plot(pvAnova, hang=-1)
pvrect(pvAnova, alpha=0.95)

## RETRIEVE MEMBERS OF SIGNIFICANT CLUSTERS.
clsig <- unlist(pvpick(pvAnova, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters)
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R") # Import tree coloring function.
dend_colored <- dendrapply(as.dendrogram(pvAnova$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
heatmap(myanova, Rowv=dend_colored, Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", RowSideColors=mycolhc)

## PLOT HEATMAP WITH HEATMAP.2() FUNCTION WHICH SCALES BETTER FOR MANY ENTRIES.
x11(height=5,width =8)
heatmap.2(myanova, Rowv=dend_colored, Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", trace="none", RowSideColors=mycolhc,margins=c(8,20))
#plot(dend_colored, horiz=T, yaxt="n")

dev.off()


## EXTRACT GENE EXPRESSION (change .txt)
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

## PRINT TABLE FOR PUBLISH
id_df_des <- data.frame(N[N$locus %in% rownames(t(NxN)),c(1,2)])
id_df_des2 <- data.frame(Locus=id_df_des$locus, ID=rep(1,188), Description=id_df_des$ID, ClusterID= paste("G", seq(1:188), sep=""), NetworkA= paste("G", seq(1:188), sep=""), NetworkB= paste("G", seq(1:188), sep=""),NetworkC= paste("G", seq(1:188), sep=""))
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\ebdbn\\")
write.table(id_df_des2, "id_anova.txt", sep="\t", quote=F, row.names=F, col.names=F)
id_def <- xtable(id_df_des2)
print(id_def, floating=T)

#########################################################################
#				SECTION -4-								#
#########################################################################
## Extract Logs ; 3 replicates per stage
## set the name and path of the locus containing file
setwd("C:/Dropbox/Workshop2013/Work/R/kegg/")
x= read.table("aatRNA.txt", sep="\t", header=F)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
NN <- genes_coc[genes_coc$locus %in% x$V1,]
cat("Dimentions of the Observed Genes Matrix :",dim(NN),"\n")

NNp <- as.matrix(NN[,4:18]) ; rownames(NNp) <- NN$locus
cat("Dimentions of the Observed Genes Matrix :",dim(NNp),"\n")

## Logs with mean per stage
library(plyr)
samples <- as.numeric(gl(5,3,15, label=seq(1:5)))
NNdf  <- data.frame(samples=samples, t(NNp))
NNm <- ddply(NNdf, .(samples), numcolwise(mean))
NNts <- t(NNm[,-1])
cat("Dimentions of the Observed Genes Matrix :",dim(NNts),"\n")

rownames(NNts) <- paste(NN[,3], NN[,2])
NNts <- NNts[order(NNts[,1]),]
windows(width=30,height=30)
heatmap.g <- heatmap(NNts, Rowv=NA, Colv=NA, col=palett, margins=c(1,40), scale="row", labCol=c("E","T","V","P","J"))


#########################################################################
#				SECTION -5-								#
#########################################################################
## SIMPLE PIECHART
## USE FOR CYTOSCAPE PROJECT
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

y <- abs(N[N$locus =="Locus554",-c(1,2)])
barplot(as.matrix(y),beside=T,col='red')



pie(y, col=rainbow(length(y), start=0.1, end=0.8), main="Pie Chart", clockwise=T)


#########################################################################
#			SECTION -1-								#
#########################################################################
require(ggplot2)
# Creates a sample data set.
## Load B2G all.set Molecular function
b2g <- read.table("C:\\Dropbox\\Workshop2013\\Work\\LaTex\\article2\\molecularFunction_all.set_cutoff1.txt", sep="\t", header=T)
colnames(b2g) <- c("id","term","score")

## Example 1
# Plots pie chart.
ggplot(b2g, aes(x = "", y = score, fill = id)) + geom_bar(width = 1) + coord_polar("y", start=pi / 3) + opts(title = "Pie Chart")
# Plots wind rose pie chart.
ggplot(b2g, aes(x = id, y = score, fill = id)) + geom_bar(width = 1) + coord_polar("y", start=pi / 3) + opts(title = "Pie Chart")

## Example 2
## Coxcomb plot = bar chart + polar coordinates
cxc <- ggplot(b2g, aes(x = score)) +
geom_bar(width = 1, colour = "black")
cxc + coord_polar()


# Windrose + doughnut plot
## both must be numeric if not skip this step
b2g$id <- cut_interval(b2g$id, length = 1)
b2g$score <- cut_number(b2g$score, 4)
## skip to here
doh <- ggplot(b2g, aes(x = id, fill = score))

# Wind rose
doh + geom_bar(width = 1) + coord_polar()
# Race track plot
doh + geom_bar(width = 0.9, position = "fill") + coord_polar(theta = "y")


## For more examples visit
http://docs.ggplot2.org/current/coord_polar.html


#########################################################################
#				SECTION -4-								#
#########################################################################
## OLD DATA
============
  Blast2GO
============
(1) write a txt file with all locus serials (extracted a priori from an anova test)
(2) load blast2go project of M edulis
(3) search and select locus serials (option: Select > select sequences by name/id)
(4) delete the rest
(5) make combined graphs
(6) export to txt data of the combined graphs
(7) create one txt file for each molecular function
(8) remove duplicate entries
(9) extract means from mev_mean.txt vs molecular function file into a sperate 3rd file
(10) draw heatmap
(11) process in Illustrator
(12) replace CUST w/ gene names

=======
  merge
=======
  ###START
  ## Merge all Locus=Names=CUST of a dataset
setwd("C:/Dropbox/Workshop2013/Work/R/anova1/")
locus <- read.table("mev_locus.txt", sep="\t", header=TRUE)
name <- read.table("mev_name.txt", sep="\t", header=FALSE) #leave it false if not you will get column duplication names
names(locus)
names(name)
ln <- merge(locus,name, by.x="locus", by.y="V1")
write.table(ln, "annotations.txt", sep="\t", quote=F) # Only once


## Variable with each dataset
## Merge the subset of a dataset. (structure is a subset of the dataset "molecular function". CUST have been extracted from
## the original dataset then merge with the file mev_name. the latter have been generated from the processing of Mev data)
means <- read.table("mev_mean_binding.txt", sep="\t", header=TRUE)
names <- ln[,2:3] # ln contains all CUST=locus=names of the dataset used to extract "means"
dtfm <- merge(means,names, by.x="name", by.y="CUST")
write.table(dtfm, "means_binding.txt", sep="\t", quote=F)
row.names(dtfm) <- dtfm$V2
dtfm # this data frame contains all means with 1 column CUST and 1 column Names (true gene names)


=============
  MEV ANOVAs
=============

## heatmaps from Anova results (MeV) 1000 permutations, Bonferoni adjusted
# Load colors
library(RColorBrewer)
?RColorBrewer
display.brewer.pal(9, "YlGnBu")
# from here
anova.d <- dtfm[!duplicated(dtfm$name),]  #remove duplicated entries
anova.d <- anova.d[order(anova.d$eggM),]  #order data from small to big in relation to egg
data <- data.matrix(anova.d)  # convert data from dataframe to matrix
data <- data[, c(2,4,6,8,10)]   #remove SDs
# make a heatmap
palett<-brewer.pal(9,"Oranges")
heatmap1 <- heatmap(data, Rowv=NA, Colv=NA, col=palett, scale="row")

====================
  Plots
====================
  ?plot
names(dtfm)
data.df <- data.frame(data) # back to a data frame
names(data.df)
data.df$name <- dtfm$name
barplot(matrix(c(data.df$eggM,data.df$veliM,data.df$juvM), nr=30, ncol=3),
        beside=TRUE, names.arg=data.df$name, col=c("red", "green", "blue"),
        width=1, xlim=c(1,100), ylim=c(-3,3), las=1, lwd=1, bty="n", cex.axis=0.6)

legend("topleft", c("Egg","Veliger","Juvenil"), pch=15, col=c("red","green","blue"), bty="n")
###END





#########################################################################
#				SECTION -2-								#
#########################################################################
## OBTAIN SIGNIFICANT CLUSTERS BY PVCLUST BOOTSTRAP ANALYSIS after a hierarchical clustering
## source : http://tinyurl.com/lqs76tf

require(pvclust)
require(gplots)
## IMPORT AND CLEANUP
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\heatmap\\")
temp <- readLines("GSE1110_series_matrix.txt")
cat(temp[-grep("^!|^\"$", temp)], file="GSE1110clean.txt", sep="\n")
mydata <- read.delim("GSE1110clean.txt", header=T, sep="\t")
rownames(mydata) <- mydata[,1]; mydata <- as.matrix(mydata[,-1])
## FILTER AND RETRIEVE HIGH INTENSITIES (50% > 100) AND HIAGH VARIABILITY (IQR > 1.5)
mydata <- mydata[apply(mydata > 100, 1, sum)/length(mydata[1,])>0.5 & apply(log2(mydata), 1, IQR) > 1.5, ]

## HIERARCHICAL CLUSTERING
mydatascale <- t(scale(t(mydata))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")

## BOOTSTRAP
pv <- pvclust(scale(t(mydata)), method.dist="correlation", method.hclust="complete", nboot=10)
plot(pv, hang=-1)
pvrect(pv, alpha=0.95)
clsig <- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters) # Retrieve members of significant clusters.
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R") # Import tree coloring function.
dend_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
heatmap(mydata, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolhc)

## PLOT HEATMAP WITH HEATMAP.2() FUNCTION WHICH SCALES BETTER FOR MANY ENTRIES.
x11(height=12)
heatmap.2(mydata, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", trace="none", RowSideColors=mycolhc)

## SORT ROWS IN DATA TABLE BY 'DEND_COLORED' AND ITS COLUMS BY 'HC'.
mydatasort <- mydata[pv$hclust$labels[pv$hclust$order], hc$labels[hc$order]]
## PLOT HEATMAP WITH BOOTSTRAP TREE IN LARGER FORMAT USING INSTEAD OF HEATMAP THE IMAGE FUNCTION.
x11(height=16, width=12)
par(mfrow=c(1,2))
plot(dend_colored, horiz=T, yaxt="n")
image(scale(t(mydatasort)), col=my.colorFct(), xaxt="n",yaxt="n")

## SAVE GRAPHICAL RESULTS TO TWO PDF FILES: 'PVCLUST.PDF' AND'HEATMAP.PDF'.
pdf("pvclust.pdf", height=21, width=6)
plot(dend_colored, horiz=T, yaxt="n")
dev.off()
pdf("heatmap.pdf", height=20, width=6)
image(scale(t(mydatasort)), col=my.colorFct(), xaxt="n",yaxt="n")
dev.off()













