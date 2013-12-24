

## Summary of all varaible levels
{
==================================
Traitement de séquences par Regexp
===================================
[L|N].*?$
[L|N].*\t.*\t
([L|N].*\t.*)(\t.*?$)
^[^C].*$
\r\n\r\n
^[C].*\t\r\n

#########################################################################
#								SECTION -1-								#
#########################################################################
========
Datasets
========
id = cust + locus (raw)	=	480 000 entrée (sondes) fromGAL file (Agilent)
id2 = cust + locus (cleaned)
exportID = locus + name + sequence =15000 séquence from Blast2GO (export)
exportID2 = locus + name (cleaned)
results = genes + stage = from limma preprocessing (annotated only for PhD)
lipid
ag
targetsISO
targetsCOC
quality
competence
performance
survival
growth

#########################################################################
#								SECTION -2-								#
#########################################################################
=======================================
Genes of interest w/ reference found in
=======================================
apoD (Moreira 2012)
MACPF proteins (calvar 2011)
P38 MAPK pathway (Irazoqui 2010 et betti 2006)
Timp et MMP genes (Tirapé 2007)
Glycoside hydrolase family 18 (GH18) includes CLPs (badariotti 2011)
STAT (canesi 2003)
BPI (zhang 2011)
Fatty acid transporter proteins (?)
Serotonin (mark fraser)
Lectin (bassem allam)
Topoisomerase (stress)
Caf (stress)
TPH Tryptophan hydroxylase (Mark fraser)
Phosphatidyl serine (http://tinyurl.com/m9qboco)
Tyrosinase (Huan 2013) (shel biogenesis)
vdg3 (venier 2006) (tissue-specific marker)

#########################################################################
#								SECTION -3-								#
#########################################################################
================================
Sampling Days post-fertilization
================================
Growth & Mortality (1, 4, 6, 8, 11, 13, 15, 18, 20, 22, 25 dpf) # Biometrie
DGE (0, 1, 12, 28, 42 dpf)	# Gene expression
Performance, Competence & Quality (0, 10, 21, 32 dpf) # profil ac gras
}







#########################################################################
#		SECTION -4-								#
#########################################################################

## Construct what=unique means, from=limma results
#merge cust + name===============================
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
id <- read.table("id2.txt", sep="\t", header=T)			## All printed probes
exportID <- read.table("Locus2ID.txt", sep="\t", header=T)			## only annotated sequences
# only annotated sequences (non annot are discarded)
cust <- merge(id,exportID)
cust <- cust[!duplicated(cust),]

#merge results + name (log ratios)
## Result or RUN files are from limma preprSocessing
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
#results <- read.table("resultsLmFitBgmin.txt", sep="\t", header=T)	## PAPER 2
results <- read.table("resultsCoc_Iso.txt", sep="\t", header=T)	## PAPER 3
results <- merge(results,cust, by.x="cust", by.y="cust")

## Consolidate duplicate rows
## calculate means of uniques probes
## either by aggregate (or with "by", or "tapply")
library(plyr)
means <- ddply(subset(results, select=-cust), .(locus), numcolwise(mean))
cat("\n","Dimensions of means before NA elimination:",dim(means),"\n")
means <- na.omit(means)
cat("\n","Dimensions of means after NA elimination:",dim(means),"\n")
means <- merge(means,exportID, by.x="locus", by.y="locus")
#description <- means[,17];means <- means[,1:16]; means <- cbind(description,means)	## PAPER 2
description <- means[,dim(means)[2]]; means <- means[,1:(dim(means)[2]-1)]; means <- cbind(description,means)	## PAPER 3
rm(description)
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\datasets\\weights\\")
write.table(means, "meansPaper3.txt", sep="\t", quote=F)
cat("\n","Dimensions of means after probe averaging:",dim(means),"\n")

#########################################################################
#								SECTION -5-								#
#########################################################################
## Extract SAM selected positive sig genes
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
siggenes_up <- read.table("samr_siggenes_up.txt", sep="\t", header=T)
locus_up <- siggenes_up[,c(3,2)]
## Extract SAM selected positive sig genes
siggenes_lo <- read.table("samr_siggenes_lo.txt", sep="\t", header=T)
locus_lo <- siggenes_lo[,c(3,2)]

#########################################################################
#								SECTION -5-								#
#########################################################################

## ANOVA MEV (not usefull)==============================================================
#anova2adj = is an anova w/ Bonferroni adjusted, and 1000 permitations
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
anova <- read.table("anovaW2.txt", sep="\t", header=F)
logs <- data.frame(gene=as.character(anova[,1]),
					locus=as.character(anova[,2]),
						Egg=anova[,3],
						Trocophore=anova[,5],
						Veliger=anova[,7],
						Pediveliger=anova[,9],
						Juvenile=anova[,11])
						#logs = genes
genes <- subset(logs, select=-Trocophore)
						genes <- as.data.frame(t(genes))

#SPecial Processing for CANOCO
d1 <- as.data.frame(subset(genes, select=-gene)); row.names(d1) <- d1$locus
d1 <- subset(d1, select=-locus); d1 <-t(d1)
		write.table(d1,"speciesW2.txt",sep="\t")

#########################################################################
#								SECTION -6-								#
#########################################################################
## Biometrie==========================================
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
growth <- read.table("Growth.txt", sep="\t", header=T)
gm <- ddply(growth, .(DPF.1), numcolwise(mean))
survival <- read.table("survival.txt", sep="\t", header=T)
sm <- ddply(survival, .(DPF.1), numcolwise(mean))

## Physiological condition
performance <- read.table("performance.txt", sep="\t", header=T)
pm <- ddply(performance, .(DPF.1), numcolwise(mean))
competence <- read.table("competence.txt", sep="\t", header=T)
cm <- ddply(competence, .(DPF.1), numcolwise(mean))
quality <- read.table("quality.txt", sep="\t", header=T)
qm <- ddply(quality, .(DPF.1), numcolwise(mean))

## Fatty acids
acids <- read.table("ag.txt", sep="\t", header=T)
fats <- read.table("lipid.txt", sep="\t", header=T)

#########################################################################
#								SECTION -7-								#
#########################################################################
## Biometrie and Fatty acid profiles (Covariables)
# Special Processing to remove additional data (so to cope with lesser data in other datasets)
## Combined dataset (0,10,21,32 dpf for Egg, veliger, pediveliger, juvenile)
stage <- c("Egg", "Veliger", "Pediveliger", "Juvenile")
# growth mean dataset where 2 stages are combined
# Case wise deletion
gm.com <- rbind(c(0,5),				# size of the egg 5µm (?)
			apply(subset(gm[gm$DPF.1 %in% c(13,15),]), 2, FUN=mean),
			apply(subset(gm[gm$DPF.1 %in% c(18,20),]), 2, FUN=mean),
			gm[gm$DPF.1==25,]); gm.com$DPF.1 <- stage

			# special processing for CANOCO
			write.table(gm.com,"growth.means.txt",sep="\t")
# survival mean combined
sm.com <- rbind(sm[sm$DPF.1==1,],
				apply(subset(sm[sm$DPF.1 %in% c(13,15),]), 2, FUN=mean),
				apply(subset(sm[sm$DPF.1 %in% c(18,20),]), 2, FUN=mean),
				sm[sm$DPF.1==25,]); sm.com$DPF.1 <- stage

pm$DPF.1 <- stage; cm$DPF.1 <- stage; qm$DPF.1 <- stage
d1 <- merge(gm.com,sm.com); d2 <- merge(d1,pm); d3 <- merge(d2,cm); physio <- merge(d3,qm)
# Fatty acids cocktail
a.com <- acids[,c(2:3,27:49)]
acids.m <- ddply(a.com, .(Traitement,Taille), numcolwise(mean))

profil <- rbind(Egg=acids.m[3,3:25],
					Veliger=acids.m[5,3:25],
					Pediveliger=acids.m[4,3:25],
					Juvenile=acids.m[2,3:25]
					)

f.com <- fats[,c(1,3,18:23)]
fats.m <- ddply(f.com, .(Traitement,Stade), numcolwise(mean))

fatty <- rbind(Egg=fats.m[4,3:8],
					Veliger=fats.m[6,3:8],
					Pediveliger=fats.m[5,3:8],
					Juvenile=fats.m[3,3:8]
					)

rm(a.com,acids,acids.m,cm,competence,cust,d1,d2,d3,anova,
exportID,f.com,fats,fats.m,gm.com,growth,id,logs,means,performance,pm
,qm,quality,results,stage,sm,sm.com,survival)





##Available datasets==================================
physio
neutral
means.W2
all.set
anova one way

## NEUTRAL**************
## Fatty acids dataset (********no trocophore, no sums*********)
{
samples <- gl(4,3,12, label=c("a.E", "b.V", "c.P", "d.J"))
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
acids <- read.table("ag.txt", sep="\t", header=T)
lipid <- read.table("lipid.txt", sep="\t", header=T)
## Percentages only of neutral lipids
neutral <- data.frame(AA=acids[c(1:6,10:12,16:18),45],
					EPA=acids[c(1:6,10:12,16:18),46],
					DHA=acids[c(1:6,10:12,16:18),47],
					SFA=acids[c(1:6,10:12,16:18),48],
					MUFA=acids[c(1:6,10:12,16:18),49],
					PUFA=acids[c(1:6,10:12,16:18),50],
					TAG=lipid[c(1:6,10:12,16:18),19],
					ST=lipid[c(1:6,10:12,16:18),21],
					AMPL=lipid[c(1:6,10:12,16:18),22],
					PL=lipid[c(1:6,10:12,16:18),23])
					}

## Profil Lipides 2011.xlsx
## PHYSIO*********************
{
library(ggplot2)
samples <- gl(4,3,12, label=c("a.E", "b.V", "c.P", "d.J"))
summary(samples)
indices <- read.table("C:\\Dropbox\\Workshop2013\\Work\\Paramètres Ecophysio\\indices.txt", sep="\t", header=T, skip=1, row.names=1)
physio <- data.frame(PI=indices[-c(7:9,13:15,19:21),7],
					OR=indices[-c(7:9,13:15,19:21),11],
					UI=indices[-c(7:9,13:15,19:21),20])

					}

## Load all GO sets
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

## Extract Logs ; 3 replicates per stage; w/ background correction | w/o between Arrays normalization | weighted (-1,2)
## EXTRACT GENE EXPRESSION
{
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
}


## Anova one way results weight 2  with background correction and betweenArrays normalization
## LOAD 5 subsets
{
setwd("C:/Dropbox/Workshop2013/Work/R/anova2/")
anova_template_EGG <- read.table("template-logs_EGG.txt", sep="\t", header=T)
anova_template_TROCO <- read.table("template-logs_TROCO.txt", sep="\t", header=T)
anova_template_VELI <- read.table("template-logs_VELI.txt", sep="\t", header=T)
anova_template_PEDI <- read.table("template-logs_PEDI.txt", sep="\t", header=T)
anova_template_JUV <- read.table("template-logs_JUV.txt", sep="\t", header=T)
lsos(pat="anova_template_.*")
}

## Transcription factors (Molecular function)
## Load Transcription factors TF.Rdata extracted from Library 1 (after molecular GO combined graph in Blast2Go)
{
setwd("C:/Dropbox/Workshop2013/Work/R/ebdbn/")
load("TF.Rdata", .GlobalEnv)
print(lsos(pat="GO.*"))

helloTF <- rbind(GO_0003700,GO_0030234,GO_0030528)
dim(helloTF)
helloTF <- helloTF[! duplicated(helloTF),]
GOx3 <- data.frame(V1=helloTF)
cat("\nDimentions of GOx3 rbind 3 GO sets before",dim(GOx3),"\n")
}

# w/ background correction to minimum| w/o between Arrays normalization | weighted (-1,2)
{					## Load SoftClustering Mfuzz package --- 12 clusters
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
load("clusters12_Lm.Rdata", .GlobalEnv)
lsos(pat="cluster.*")
}

## Logs of 14,876 transcript after pre-processing
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
N <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)	## W/n+Weight (PAPER 2)
N <- read.table("meansW2.txt", sep="\t", header=T)		## Bg+W/n+Betw+Weight

## Library 1 for subset A, B, C w background correction to minimum | w/o between Arrays normalization | weighted (-1,2)
{
setwd("C:/Dropbox/Workshop2013/Work/R/sam/")
all.set <- read.table("allSet_Lm.txt", sep="\t", header=T)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/weights/")
genes_coc <- read.table("meansLmFitBgmin.txt", sep="\t", header=T)
genes_SIG <- genes_coc[genes_coc$locus %in% all.set$x,]
rownames(genes_SIG) <- genes_SIG$locus
species <- as.matrix(t(genes_SIG[,-c(1,2)]))
dim(species)	#15x3633
}

## Anova one way results w/ background correction to minimum | w/o between Arrays normalization | weighted (-1,2)
## LOAD 5 subsets (from 5 different Templates)
{
setwd("C:/Dropbox/Workshop2013/Work/R/ANOVA/")
anova_template_EGG <- read.table("meansLmFitBgmin_EGG.txt", sep="\t", header=T)
anova_template_VELI <- read.table("meansLmFitBgmin_VEL.txt", sep="\t", header=T)
anova_template_JUV <- read.table("meansLmFitBgmin_JUV.txt", sep="\t", header=T)
lsos(pat="anova_template_.*")
}
{
setwd("C:/Dropbox/Workshop2013/Work/R/ANOVA/")
load("anova.Rdata", .GlobalEnv)
lsos(pat="nova")
}


## Transcription factors, kinases, phosphatases, chromatin remodelers
## FOR MORE INFO HEAD TO ENRICHMENT
{
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
load("regulators.Rdata", .GlobalEnv)
lsos(pat="GO2gene_.*|.*_matrix|.*2fold|TF16")
}

## SUBSET E
{
setwd("C:/Dropbox/Workshop2013/Work/R/enrichment/")
load("subsetE2.Rdata", .GlobalEnv)
lsos(pat="GO2gene_.*|.*_matrix|.*setE")
## GROWTH (G), DEVELOPMENT (D), LOCOMOTION (L), RHYTHMIC (R)
## Neural(Ne), Tissue (Ti), Muscle (Mu)
## All_matrix & All_matrix2 : combined genes of every search
## logsSubsetE & logsSubsetE2 : logs of genes before 2fold filtering
## sig_All_matrix : combined genes found in Library1
}

## NETWORKS EBDBN
setwd("C:/Dropbox/Workshop2013/Work/R/ebdbn/")
load("Networks2.3.Rdata", .GlobalEnv)
lsos(pat="net.*E.|.*D1|.*B1|networkF.")




