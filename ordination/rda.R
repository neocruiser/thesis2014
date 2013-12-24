## Tutorials
# RDA - vegantutor
# RDA - dissimilarities
# RDA - vegan package
# RDA - Legendre & Legendre 1998
# RDA - intro-vegan

#########################################################################
#								SECTION -1-								#
#########################################################################
===================
Redundancy Analysis FAT and ENV (dummy variable) 	-- CONSTRAINED --
===================
## RDA & rdaTest
## Legendre

## Work with 3 replicates COCktail
## Fatty acids dataset (no trocophore, no sums)
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
acids <- read.table("ag.txt", sep="\t", header=T)
fat <- acids[c(1:6,10:12,16:18),27:46]
# Centring by col
species <- NULL
for(j in 1:12) species <- rbind(species, (fat[j,]-(apply(fat,2,mean))))
species
}

## Constraining by a categorical variable
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
env <- read.table("env.txt", sep="\t", header=T)
env.ns <- env[-c(4:6),-c(2)]
}

## Standardize to mean 0 and variance 1
## ************Method (2)*********** (method 1 cf section 2)
{
env <- NULL
for(j in 1:12) env <- rbind(env, (env.ns[j,]-(cov.wt(env.ns)$center))/(apply(env.ns,2,sd)))
env
env <- as.matrix(env)	
}		


## RUN an RDA (linearly constrained PCA)	
library(vegan)			
result <- rda(species, env)
summary(result)
# Extract the species and site scores on the ordination axes 
#    choice=c(1,2) indicates the first two axes
#    scaling=1 (site), 2 (species), or 3 (both) are scales by eigenvalues
const.sites <- sqrt((nrow(species)-1)*sum(eigenvals(result)))
scores(result, choice=c(1,2), scaling=1, const=c(1,const.sites))

plot(result, dis=c("cn","sp", "wa"), type="n", xaxt="n", yaxt="n")
axis(1,col.axis="forestgreen", col.ticks="forestgreen")
axis(2,col.axis="forestgreen", col.ticks="forestgreen")
axis(3,col.axis="brown", col.ticks="brown")
axis(4,col.axis="brown", col.ticks="brown")
text(result, dis="cn", col="chocolate")
points(result, display="sp",pch=21, bg="forestgreen", cex=0.8)
text(result, display="wa",font=1, cex=0.7)
#plot(envfit(result, species, display="lc"), add = TRUE, col="red")

## Legendre Package
library(rdaTest)
result2 <- rdaTest(species, env, test.F=FALSE, nperm=999)
summary(result2)

# distance biplot	(binary values)
plot(result2, graph.type="Z",scaling=1)
# correlation biplot 			(quantitativ values)
plot(result2, graph.type="Z",scaling=2)



#########################################################################
#								SECTION -2-								#
#########################################################################
## RDA of GENES and FAT
## RDA & rdaTest
## Legendre
## Work with 3 replicates Cocktail microarrays
## Genes
## all.set for subset A, B, C
{
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
all.set <- read.table("all.set.txt", sep="\t", header=F)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansW2.txt", sep="\t", header=T)
genes_SIG <- genes_coc[genes_coc$locus %in% all.set$V2,]
rownames(genes_SIG) <- genes_SIG$locus
logs <- as.matrix(t(genes_SIG[,4:18]))
dim(logs)	#15x1444
species <- logs[-c(4:6),]
}
## LOAD neutral and polar lipids
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
acids <- read.table("ag.txt", sep="\t", header=T)
lipid <- read.table("lipid.txt", sep="\t", header=T)
## Percentages only of neutral lipids
dataset <- data.frame(AA=acids[c(1:6,10:12,16:18),45],
					EPA=acids[c(1:6,10:12,16:18),46],
					DHA=acids[c(1:6,10:12,16:18),47],
					SFA=acids[c(1:6,10:12,16:18),48],
					MUFA=acids[c(1:6,10:12,16:18),49],
					PUFA=acids[c(1:6,10:12,16:18),50],
					TAG=lipid[c(1:6,10:12,16:18),19],
					ST=lipid[c(1:6,10:12,16:18),21],
					AMPL=lipid[c(1:6,10:12,16:18),22],
					PL=lipid[c(1:6,10:12,16:18),23])
## Load Indices
indices <- read.table("C:\\Dropbox\\Workshop2013\\Work\\Paramètres Ecophysio\\indices.txt", sep="\t", header=T, skip=1, row.names=1)
physio <- data.frame(PI=indices[-c(7:9,13:15,19:21),7], 
					OR=indices[-c(7:9,13:15,19:21),11],
					UI=indices[-c(7:9,13:15,19:21),20])
					}
## Standardize to mean 0 and variance 1
## ************Method (1)***********
{
require(vegan)
species <- decostand(species, method="range")
env <- decostand(dataset, method="standardize")		# dataset or physio
}

## RUN an RDA (linearly constrained PCA)		
result <- rda(species, dataset)
{
plot(result, dis=c("cn","sp"), type="n", xaxt="n", yaxt="n")
axis(1,col.axis="forestgreen", col.ticks="forestgreen")
axis(2,col.axis="forestgreen", col.ticks="forestgreen")
axis(3,col.axis="brown", col.ticks="brown")
axis(4,col.axis="brown", col.ticks="brown")
text(result, dis="cn", col="chocolate")
points(result, display="sp",pch=21, bg="forestgreen", cex=0.8)
}

## Legendre Package
{
lipids <- dataset[,c(5,6,7,9)]
env <- decostand(lipids, method="standardize")
library(rdaTest)
result2 <- rdaTest(species, lipids, test.F=T, nperm=999)
summary(result2)
result2$Rsquare
}
?rdaTest
# distance biplot		(binary values)
plot(result2, graph.type="Z",scaling=1, label.spe=F)
# correlation biplot			(quantitativ values)
plot(result2, graph.type="Z",scaling=2, plot.spe=F)
points(result2, plot.spe=T, label.spe=F)


#########################################################################
#								SECTION -3-								#
#########################################################################
## PAPER 2
set.seed(1593058)
require(vegan)

{					## Load L.Rdata
setwd("C:/Dropbox/Workshop2013/Work/R/ordination/")
load("rda.Rdata", .GlobalEnv)
lsos(pat="")
}
## Library 1 (2nd paper)
## Allspecies=w/ troco || use species=w/o troco only for lipids =====IMPORTANT=====
{
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
all.set <- read.table("allSet_Lm.txt", sep="\t", header=T)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
genes_SIG <- genes_coc[genes_coc$locus %in% all.set$x,]
rownames(genes_SIG) <- genes_SIG$locus
Allspecies <- as.matrix(t(genes_SIG[,-c(1,2)]))
species <- Allspecies[-c(4:6),]			## w/o trocophores
dim(Allspecies)	#15x3633						## w/ trocophores
}
## LOAD neutral and polar lipids
{
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
acids <- read.table("ag.txt", sep="\t", header=T)
lipid <- read.table("lipid.txt", sep="\t", header=T)
## Percentages only of neutral lipids
dataset <- data.frame(AA=acids[c(1:6,10:12,16:18),45],
					EPA=acids[c(1:6,10:12,16:18),46],
					DHA=acids[c(1:6,10:12,16:18),47],
					SFA=acids[c(1:6,10:12,16:18),48],
					MUFA=acids[c(1:6,10:12,16:18),49],
					PUFA=acids[c(1:6,10:12,16:18),50],
					TAG=lipid[c(1:6,10:12,16:18),19],
					ST=lipid[c(1:6,10:12,16:18),21],
					AMPL=lipid[c(1:6,10:12,16:18),22],
					PL=lipid[c(1:6,10:12,16:18),23])
## Load Indices
indices <- read.table("C:\\Dropbox\\Workshop2013\\Work\\Paramètres Ecophysio\\indices.txt", sep="\t", header=T, skip=1, row.names=1)
physio <- data.frame(PI=indices[-c(7:9,13:15,19:21),7], 
					OR=indices[-c(7:9,13:15,19:21),11],
					UI=indices[-c(7:9,13:15,19:21),20])
					}
## Create binary matrix
## AllbinaryStage=w/ troco || binaryStage=w/o troco for lipids ===== IPORTANT=====
## Binary for DEvelopmental stages (w/ trocophore)
{
stage <- as.data.frame(model.matrix(~0+factor(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3)))))
colnames(stage)=c("Egg","Trocophore","Veliger","Pediveliger","Juvenile")
stage
}
## Binary for DEvelopmental stages (w/o trocophore)
{
stage <- as.data.frame(model.matrix(~0+factor(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3)))))
colnames(stage)=c("Egg","Veliger","Pediveliger","Juvenile")
## Create hybrid matrix of lipid and binary (only 12 rows, no trocophores)

hybrid <- cbind(dataset,stage)
}
## Standardize (optional)
{
spec <- decostand(species, method="range")		# method=log combined w log(lipids) produces an interesting RDA w/o troco
Allspec <- decostand(Allspecies, method="range")		# method=log w/ troco w/o (lipids)
lipids <- decostand(dataset, method="range")	#method=standardize produces same results; log produces a slightly different results (but same conclusions)
phys <- decostand(physio, method="standardize")
env <- decostand(as.matrix(cbind(dataset,physio)), method="range")
hybrid <- decostand(hybrid,method="range") ## w/o troco w/lipids
}
## Load all GO sets (From Fishers enrichment analysis) -- optional
{
GO1 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0043232.txt", sep="\t", header=F)
GO2 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0006412.txt", sep="\t", header=F)
GO3 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0010467.txt", sep="\t", header=F)
GO4 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0030529.txt", sep="\t", header=F)
GO5 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0005694.txt", sep="\t", header=F)
GO6 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0005198.txt", sep="\t", header=F)
GO7 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0032991.txt", sep="\t", header=F)
GO8 <- read.table("C:\\Dropbox\\Workshop2013\\Work\\R\\ordination\\GO_0044424.txt", sep="\t", header=F)
lsos(pat="GO.")
}
## Load all clusters from (last soft clustering for paper 2 w all.set)
{					## Load SoftClustering Mfuzz package --- 12 clusters
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
load("clusters12_Lm.Rdata", .GlobalEnv)
lsos(pat="cluster.*")
}
## Transcription factors (Molecular function)
## Load Transcription factors TF.Rdata extracted from Library 1 (after molecular GO combined graph in Blast2Go)
{
## Transcription factors, kinases, phosphatases, chromatin remodelers
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
load("regulators.Rdata", .GlobalEnv)
lsos(pat="GO2gene_.*|.*_matrix|.*2fold|TF16")

#	GOx3 <- data.frame(V1=TF2fold$locus)
GOx3 <- data.frame(V1=TF16$Locus)
cat("\nDimentions of GOx3 rbind 3 GO sets before",dim(GOx3),"\n")
}



============
## RDA + Lipids
============
result <- rda(spec ~ PUFA+TAG+MUFA+AMPL,lipids)
{				## diagnostics
vif.cca(result)
result$CCA$tot.chi/result$tot.chi			## Total variance
RsquareAdj(result)
anova(result, step=2000)				# by= term; axis; margin (just in ~formula)
plot(procrustes(rda(spec), extended))
## Test of significance
envfit(result~ TAG+AMPL+PUFA, lipids, perm=999, dis="lc", scaling=2)
}
## PLOTTING
a=c(1:16)			## Choose which Transription factors to DRAW
{				## Empty plot
x <- scores(result)$species
pl <- plot(result,dis=c("cn","sp"), xaxt="n", yaxt="n", scaling=2, type="n")
points(x[rownames(x) %in% GOx3[a,1],], pch=20, col="grey", cex=1)		## Transcription factors (run seperatly then run "identify"...)
axis(1,col.axis="forestgreen", col.ticks="forestgreen")
axis(2,col.axis="forestgreen", col.ticks="forestgreen",las=1)
axis(3,col.axis="brown", col.ticks="brown")
axis(4,col.axis="brown", col.ticks="brown",las=1)

points(result, dis="sp", pch=1, col="grey")
}
{				## Draw points
## Extract species scores (for additional manipulation see section 1)
x <- scores(result)$species
## create color palette
y <- rainbow(30)
points(x[rownames(x) %in% cluster1[,1],], pch=20, col=y[1], , cex=1.1)
points(x[rownames(x) %in% cluster2[,1],], pch=20, col=y[22], , cex=1.1)
points(x[rownames(x) %in% cluster3[,1],], pch=20, col="black", cex=1.1)
points(x[rownames(x) %in% cluster4[,1],], pch=20, col=y[13], cex=1.1)
points(x[rownames(x) %in% cluster5[,1],], pch=20, col=y[18])
points(x[rownames(x) %in% cluster6[,1],], pch=20, col=y[24])
points(x[rownames(x) %in% cluster7[,1],], pch=20, col=y[7])
points(x[rownames(x) %in% cluster8[,1],], pch=20, col=y[16], cex=1)
points(x[rownames(x) %in% cluster9[,1],], pch=20, col=y[25])
points(x[rownames(x) %in% cluster10[,1],], pch=20, col=y[28])
points(x[rownames(x) %in% cluster11[,1],], pch=20, col=y[4] , cex=1)
points(x[rownames(x) %in% cluster12[,1],], pch=20, col=y[6], cex=1.1)

points(x[rownames(x) %in% GOx3[,1],], pch=1, col="black", cex=1)					## show all transcription factors selected in GOx3
#segments(0, 0, x[rownames(x) %in% GOx3[a,1],1], x[rownames(x) %in% GOx3[a,1],2], col="forestgreen", lwd="0.7")			## plot lines to interesting TF
text(result, dis="cn", col="chocolate", font=4)
}
## Interactive identity function
## Hand picked annotation on the graph
identify(pl, "sp", labels=colnames(spec), cex=0.55, col="black")		

## Model Building with "step" + AIC
## Reducing the number of constraints (VeganTutor.pdf)
{
cat("\n\nIf constrained inertia in permutation test is lower than the observed inertia: Constraints are then significant")
## Start with unconstrained (option 1)
extended <- rda(spec ~. , lipids)			# Dont Manipulate
reduce0 <- rda(spec ~1, lipids)				# ~1 for unconstrained
reduce <- step(reduce0, scope=formula(extended), test="perm")
reduce
reduce$anova
vif.cca(reduce)
## start with constraint model (option 2)
reduce.all <- step(extended, scope=list(lower=formula(reduce0), upper=formula(extended)),trace=0)
reduce.all
reduce.all$anova
vif.cca(reduce.all)
}

## 3D plotting
{
library(rgl)
library(scatterplot3d)
ordirgl(result, size="2")
ordirgl(result, display="species", type="p", col="yellow")
## 3D plotting (type 2)
ol <- ordiplot3d(result, scaling = 3, angle=15, type="n")
points(ol, "points", pch=16, col="red", cex = 0.7)
}

## OPTIONAL (label only specific set, and draw points for everything else)
orditorp(result, disp="sp",col="violetred2", pcol="grey", pch="")
orditorp(result, disp="sp",priority=GOx3[,1], col="forestgreen", pcol="violetred2", pch=20, cex=0.65)
orditorp(result, disp="sp", labels=c(rownames(x[rownames(x) %in% GOx3[,1],])),select=GOx3[,1], pch=20, pcol="blue", cex=0.9)

## Save a copy of the memory
setwd("C:/Dropbox/Workshop2013/Work/R/ordination/")
lsos(pat="")
save(list=ls(pattern=""),file="rda.Rdata")

============
## RDA + Binary
============
## Allspec= w/ troco w/ lipids
result <- rda(Allspec~Egg+Veliger+Pediveliger+Juvenile,stage)
## PLOTTING
a=c(1:16)			## Choose which Transription factors to DRAW
{				## Empty plot
x <- scores(result)$species
pl <- plot(result,dis=c("cn","sp"), xaxt="n", yaxt="n", scaling=2, type="n")
points(x[rownames(x) %in% GOx3[a,1],], pch=20, col="grey", cex=1)		## Transcription factors (run then stop then run "identify"...)
axis(1,col.axis="forestgreen", col.ticks="forestgreen")
axis(2,col.axis="forestgreen", col.ticks="forestgreen",las=1)
axis(3,col.axis="brown", col.ticks="brown")
axis(4,col.axis="brown", col.ticks="brown",las=1)

points(result, dis="sp", pch=1, col="grey",cex=0.5)
}
{				## Draw points
## Extract species scores (for additional manipulation see section 1)
x <- scores(result)$species
## create color palette
#y <- rainbow(30) 		## then use col=y[1]
points(x[rownames(x) %in% cluster1[,1],], pch=20, col="#028E9B", , cex=1.1)
points(x[rownames(x) %in% cluster2[,1],], pch=20, col="#133CAC", , cex=1.1)
points(x[rownames(x) %in% cluster3[,1],], pch=20, col="#25D500", cex=1.1)
points(x[rownames(x) %in% cluster4[,1],], pch=20, col="#FF6A00", cex=1.1)
points(x[rownames(x) %in% cluster5[,1],], pch=20, col="#7C07A9",, cex=0.92)
points(x[rownames(x) %in% cluster6[,1],], pch=20, col="#D30068")
points(x[rownames(x) %in% cluster7[,1],], pch=20, col="black")
points(x[rownames(x) %in% cluster8[,1],], pch=20, col="#F0FC00", cex=0.92)
points(x[rownames(x) %in% cluster9[,1],], pch=20, col="#FF4900")
points(x[rownames(x) %in% cluster10[,1],], pch=20, col="#04B8E2", cex=0.92)
points(x[rownames(x) %in% cluster11[,1],], pch=20, col="#E40045" , cex=1)
points(x[rownames(x) %in% cluster12[,1],], pch=20, col="#00AF64", cex=1.1)

points(x[rownames(x) %in% GOx3[,1],], pch=1, col="black", cex=1)					## show all transcription factors selected in GOx3
#segments(0, 0, x[rownames(x) %in% GOx3[a,1],1], x[rownames(x) %in% GOx3[a,1],2], col="forestgreen", lwd="0.7")			## plot lines to interesting TF
text(result, dis="cn", col="chocolate", font=4)
}
## Interactive identity function
## Hand picked annotation on the graph
identify(pl, "sp", labels=colnames(spec), cex=0.55, col="black")		
{				## diagnostics
vif.cca(result)
result$CCA$tot.chi/result$tot.chi			## Total variance
RsquareAdj(result)
anova(result, step=2000)				# by= term; axis; margin (just in ~formula)
	#plot(procrustes(rda(spec), extended))
## Test of significance
envfit(result~Egg+Trocophore+Veliger+Pediveliger+Juvenile, stage, perm=999, dis="lc", scaling=2)
}
## Model Building with "step" + AIC
## Reducing the number of constraints (VeganTutor.pdf)
{
cat("\n\nIf constrained inertia in permutation test is lower than the observed inertia: Constraints are then significant")
## Start with unconstrained (option 1)
extended <- rda(Allspec ~. , stage)			# Dont Manipulate
reduce0 <- rda(Allspec ~1, stage)				# ~1 for unconstrained
reduce <- step(reduce0, scope=formula(extended), test="perm")
reduce
reduce$anova
vif.cca(reduce)
## start with constraint model (option 2)
reduce.all <- step(extended, scope=list(lower=formula(reduce0), upper=formula(extended)),trace=0)
reduce.all
reduce.all$anova
vif.cca(reduce.all)
}











