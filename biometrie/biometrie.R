## IMPORTANT***
## Load the summarySE function from 01funcs

========================
  Physiological studies
========================

## load datasets
setwd("C:/Dropbox/Workshop2013/Work/Paramètres Ecophysio/")
load("biometrics.Rdata", .GlobalEnv)
lsos(pat="")

## DATASET 1
## source: C:\Dropbox\Workshop2013\Work\Paramètres Ecophysio\Profil Lipides 2011.xlsx (feuillet: organized)
## Indices (Peroxidation index & ratio n3/n6)
comp.indices <- read.table("clipboard", sep="\t", header=T)
indices <- comp.indices[-c(1:3),]
indices

## DATASET 2
## NEUTRAL (AA, EP, DHA + Iatro lipids)
## Fatty acids dataset (no trocophore, no EGG)
fattyacids <- read.table("clipboard",sep="\t", header=T)
fattyacids
neutral <- fattyacids[-c(1:9),c(3:4,6:15)]
neutral

## in the case of repeated measurements on the same group (e.g., of animals, individuals, cultures, or reactions),
## CIs or SE bars are irrelevant to comparisons within the same group
## Read biometrie file
setwd("C:/Dropbox/Workshop2013/Work/R/Biometrie/")
biometrie <- read.table("Biometrie.txt", sep="\t", header=T)
names(biometrie)
library(ggplot2)
windows(record=T)

## GROWTH for both diets COC and TISO
sgrowth <- biometrie[, c(4,6,8)]
names(sgrowth)
sgrowth <- sgrowth[complete.cases(sgrowth),]   # Remove NA rows
		# write.table(sgrowth, "sgrowth.txt", sep="\t")
DF <- summarySE(sgrowth, measurevar="Growth", groupvars=c("DPF.1","Treatment"))
DF
                                        # Confidence interval of the mean of Mortality
ggplot(DF, aes(x=DPF.1, y=Growth)) +
  geom_errorbar(aes(ymin=Growth-se, ymax=Growth+se), width=.4) +
  geom_line(aes(group=Treatment)) +
  geom_point(size=4.4, fill="black",aes(shape=Treatment)) +
    theme_bw()

## Mortality COC and TISO
names(biometrie)
survival <- biometrie[, c(4,6,9)]
names(survival)
?complete.cases
survival <- survival[complete.cases(survival),]   # Remove NA rows
		# write.table(survival, "survival.txt", sep="\t")
DF <- summarySE(survival, measurevar="Mortality", groupvars=c("DPF.1","Treatment"))
# Confidence interval of the mean of Mortality
ggplot(DF, aes(x=DPF.1, y=Mortality)) +
  geom_errorbar(aes(ymin=Mortality-se, ymax=Mortality+se), width=.4) +
  geom_line(aes(group=Treatment)) +
  geom_point(size=4.4, fill="black",aes(shape=Treatment)) +
    theme_bw()


====================
    Summary Stats
====================
## PLOT Neutral(AA, EPA, DHA ...., TAG, PL, ST ....) and indices (OR, PI)

## Plot INDICES
samples <- gl(3,3,18, label=c("1Veliger", "2Pediveliger", "3Juvenile"))
diet <- gl(2,9,18, label=c("COC","Tiso"))
DF <- summarySE(indices[,-c(1:2)], measurevar="PI", groupvars=c("samples","diet"))
DF
## Plot SEM (standard error of the mean)
## source : http://tinyurl.com/lo53qls
require(ggplot2)
pd <- position_dodge(.1)
ggplot(DF, aes(x=samples, y=PI, colour='black', group=diet))+
    geom_errorbar(aes(ymin=PI-se, ymax=PI+se), colour='black', width=.1, position=pd)+
    geom_line(position=pd, colour="black")+
    geom_point(position=pd, size=5,aes(shape=diet), colour="black")+
    xlab("Stages of development (DPF)")+
    ylab("PI")+
    scale_colour_hue(name="Treatment",
                     breaks=c("COC", "Tiso"),
                     labels=c('Cocktail', 'Tiso'),
                     l=40)+
    ggtitle("The effect of a fatty acid deficit\ntreatment on physiological indices")+
    scale_y_continuous(limits=c(0, max(DF$PI + DF$se)),				## set y range
#                             breaks=.25:25*0.2)+				## for OR ... set tick every .25
    			breaks=0:15*25)+				## for PI ... set tick every 25
    theme_bw()+
    theme(legend.justification=c(1,0), legend.position=c(1,0))

## OR
DF <- summarySE(indices[,-c(1:2)], measurevar="OR", groupvars=c("samples","diet"))
pd <- position_dodge(.1)
ggplot(DF, aes(x=samples, y=OR, colour='black', group=diet))+
    geom_errorbar(aes(ymin=OR-se, ymax=OR+se), colour='black', width=.1, position=pd)+
    geom_line(position=pd, colour="black")+
    geom_point(position=pd, size=5,aes(shape=diet), colour="black")+
    xlab("Stages of development (DPF)")+
    ylab("OR")+
    scale_colour_hue(name="Treatment",
                     breaks=c("COC", "Tiso"),
                     labels=c('Cocktail', 'Tiso'),
                     l=40)+
    ggtitle("The effect of a fatty acid deficit\ntreatment on physiological indices")+
    scale_y_continuous(limits=c(0, max(DF$OR + DF$se)),				## set y range
                             breaks=.25:25*0.2)+				## for OR ... set tick every .25
    theme_bw()+
    theme(legend.justification=c(1,0), legend.position=c(1,0))


## Plot neutral elements
## Prepare dataset
samples <- gl(3,3,18, label=c("1Veliger", "2Pediveliger", "3Juvenile"))
diet <- gl(2,9,18, label=c("COC","Tiso"))
DFnall <- NULL
for(i in 1:(dim(neutral)[2]-2)){
    DF <- summarySE(neutral[,-c(1:2)], measurevar=names(neutral[,-c(1:2)])[i], groupvars=c("samples","diet"))
    names(DF)[names(DF)== names(neutral[,-c(1:2)])[i]] <- "effect"			## change names of the varaible column
    DF <- cbind(DF, variable=names(neutral[,-c(1:2)])[i])
    DFnall <- rbind(DFnall, DF)
}
DFnall
## the loop creates a single dataframe of the summary of all variables in neutral

## plot
## just change the DFnall rows from DFnall[1:18,] to DFnall[19:36,] to DFnall[37:60,]
pd <- position_dodge(.1)
ggplot(DFnall[37:60,], aes(x=samples, y=effect,  group=interaction(diet,variable)))+
    geom_errorbar(aes(ymin=effect-se, ymax=effect+se), colour='black', width=.3, position=pd)+
    geom_line(position=pd, aes(colour=diet))+
    geom_point(position=pd, size=4.4,aes(shape=variable), colour="black")+
    xlab("Stages of development (DPF)")+
    ylab("effect")+
    scale_colour_hue(name="Treatment",
                     breaks=c("COC", "Tiso"),
                     labels=c('Cocktail', 'Tiso'),
                     l=40)+
    scale_y_continuous(limits=c(0, max(DFnall[37:60,4] + DFnall[37:60,6])),				## set y range
    			breaks=0:35*2)+				## set tick every 25
    ggtitle("The effect of a fatty acid deficit\ntreatment on physiological indices")+
    theme_bw()+
    theme(legend.justification=c(1,0), legend.position=c(1,0))


## plot using lattice
require(lattice)
xyplot(effect ~ samples | diet + variable, data=DFnall, type="o", col="black",pch=20, cex=1)

## plot with scale activated
## Plot used for PAPER 3
require(ggplot2)
pd <- position_dodge(.1)
ggplot(DFnall, aes(x=samples, y=variable, group=diet)) +
    geom_point(aes(size=effect, colour=diet), shape=22, position=pd) +
    #geom_point(aes(size=effect-se), colour='black', position=pd, alpha=I(.2)) +
    scale_size_area(max_size=10) +
    scale_colour_hue(name="Treatment",
                     breaks=c("COC", "Tiso"),
                     labels=c('Cocktail', 'Tiso'),
                     l=40)+
    ggtitle("The effect of a fatty acid deficit\ntreatment on physiological lipids")+
    theme_bw()+
    theme(legend.justification=c(1,0), legend.position=c(1,0))+
    xlab("Stages of development (DPF)")+
    ylab("Fatty acids")

## Plot neutral EGG stage only
## Prepare dataset
## Fatty acids dataset (EGG)
fattyacids <- read.table("clipboard",sep="\t", header=T)
fattyacids
neutral <- fattyacids[-c(1:6),c(3:4,6:15)]
neutral

samples <- gl(1,3,3, label=c("Egg"))
diet <- gl(1,3,3, label=c("none"))
DFnall <- NULL
for(i in 1:(dim(neutral)[2]-2)){
    DF <- summarySE(neutral[,-c(1:2)], measurevar=names(neutral[,-c(1:2)])[i], groupvars=c("samples","diet"))
    names(DF)[names(DF)== names(neutral[,-c(1:2)])[i]] <- "effect"			## change names of the varaible column
    DF <- cbind(DF, variable=names(neutral[,-c(1:2)])[i])
    DFnall <- rbind(DFnall, DF)
}
DFnall
## the loop creates a single dataframe of the summary of all variables in neutral

require(ggplot2)
pd <- position_dodge(.1)
ggplot(DFnall, aes(x=samples, y=variable, group=diet)) +
    geom_point(aes(size=effect, colour=diet), shape=22, position=pd) +
    #geom_point(aes(size=effect-se), colour='black', position=pd, alpha=I(.2)) +
    scale_size_area(max_size=10) +
    scale_colour_hue(name="Treatment",
                     breaks=c("COC", "Tiso"),
                     labels=c('Cocktail', 'Tiso'),
                     l=40)+
    ggtitle("The effect of a fatty acid deficit\ntreatment on physiological lipids")+
    theme_bw()+
    theme(legend.justification=c(1,0), legend.position=c(1,0))+
    xlab("Stages of development (DPF)")+
    ylab("Fatty acids")


## save grouped datasets
setwd("C:/Dropbox/Workshop2013/Work/Paramètres Ecophysio/")
lsos(pat="indices|fatty.*|neutral|biometrics")
save(list=ls(pattern="indices|fatty.*|neutral|biometrics"),file="biometrics.Rdata")

====================
    ## + Diagnostics
====================

## Diagnostics (Anova, Multiple comparison)
## one way ANOVA (2 or more conditions)
## START (run)
cat("\n R2 represents the % of variance in the response variable accounted for or explained by the variability in the predictor variable\n\n")

## Anova & Tukey HSD multicomparison test
diet <- gl(2,9,18, label=c("COC","Tiso"))
samples <- gl(3,3,18, label=c("1Veliger", "2Pediveliger", "3Juvenile"))

## All samples + Both treatments
anova.neutral.all <- apply(neutral[,-c(1:2)], 2, function(x) {
list(
fit <- aov(x ~ diet + samples + diet*samples),
p <- drop1(fit, ~., test="F"),
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),
TukeyHSD(fit, ordered=T)
)})

## All samples + one treatment
## COC treatment
apply(neutral[1:9,-c(1:2)], 2, function(x) {
list(
fit <- aov(x ~ samples[1:9]),
p <- drop1(fit, ~., test="F"),
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),
TukeyHSD(fit, ordered=T)
)})

## Tiso treatment
apply(neutral[10:18,-c(1:2)], 2, function(x) {
list(
fit <- aov(x ~ samples[1:9]),
p <- drop1(fit, ~., test="F"),
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),
TukeyHSD(fit, ordered=T)
)})

## All samples + Both treatments
apply(indices, 2, function(x) {
list(
fit <- aov(x ~ samples + diet + samples*diet),
p <- drop1(fit, ~., test="F"),
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),
TukeyHSD(fit, ordered=T)
)})

## All samples + one treatment
## COC treatment
apply(indices[1:9], 2, function(x) {
list(
fit <- aov(x ~ samples[1:9,]),
p <- drop1(fit, ~., test="F"),
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),
TukeyHSD(fit, ordered=T)
)})

## All samples + one treatment
## Tiso treatment
apply(indices[10:18], 2, function(x) {
list(
fit <- aov(x ~ samples[1:9,]),
p <- drop1(fit, ~., test="F"),
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),
TukeyHSD(fit, ordered=T)
)})


## Bonferroni correction (multiple comparison) -- It can be used only if the predictor is a 2 levels group (which is only the case for multiple comparison between stages of the same diet treatment= not useful for paper3)
Neutral_nonMultipleComparison <- apply(neutral[-c(1:2)], 2, function(x) pairwise.t.test(x, samples, p.adj= "none"))
Neutral_bonferroniMultipleComparison <- apply(neutral, 2, function(x) pairwise.t.test(x, a, p.adj= "bonf"))
Physio_bonferroniMultipleComparison <- apply(physio, 2, function(x) pairwise.t.test(x, a, p.adj= "bonf"))

	## Dunnett test (more conservative)
	require(multcomp)
	?glht
	fit.dunnett <- glht(fit,lincft=mcp(method="Dunnett"))
	confint(fit.dunnett)
	plot(fit.dunnett)

## MEan + SD
tapply(y,x,mean);tapply(y,x,sd)

## Plot anova summary
layout(matrix(c(1:4),2,2))
x11(); plot(fit)
## One factor plotting
library(gplots)
## Plot 2
x(11); plotmeans(y ~ x, connect=list(1:2,3:4),bars=T,lwd=2, xlab="Development stages", ylab="")
## Plot 3 (boxplot)
x11(); boxplot.n(y ~ x, shrink=0.8, notch=F, boxwex=0.25, las=1)
## END (run)e

## Distribution (Shapiro-Wilk | skewness | ...)e
library(e1071)
neutral_skewed <- apply(neutral, 2, function(x) skewness(x))
physio_skewed <- apply(physio, 2, function(x) skewness(x))
neutral_kurtosis <- apply(neutral, 2, function(x) kurtosis(x))
physio_kurtosis <- apply(physio, 2, function(x) kurtosis(x))

## Shapiro wilk (W). Test for normality
cat("\n\n", "Test for normality. \nIf p-value is significant the null hypothesis is rejected, and tentatively conclude that X is not normally distributed. \nX displays a statistically significant departure from normality\n")
neutral_normality <- apply(neutral, 2, function(x) shapiro.test(x))
physio_normality <- apply(physio, 2, function(x) shapiro.test(x))

## STEM and leaf plot
neutral_stem <- apply(neutral, 2, function(x) stem(x))
physio_stem <- apply(physio, 2, function(x) stem(x))

## Kolmogorov-Smirnov test
neutral_normality_kolSmir <- apply(neutral, 2, function(x) ks.test(x, "pnorm", mean=mean(x), sd=sqrt(var(x))))
physio_normality_kolSmir <- apply(physio, 2, function(x) ks.test(x, "pnorm", mean=mean(x), sd=sqrt(var(x))))

## HISTOGGRAM and density lines
par(mfrow=c(3,4))
for(i in 1:ncol(neutral)){
hist(neutral[,i], main=paste("Histogram of",names(neutral[i])), xlab=names(neutral[i]), las=1)
lines(density(neutral[,i], bw=0.2))
rug(neutral[,i])
}

par(mfrow=c(1,3))
for(i in 1:ncol(physio)){
hist(physio[,i], main=paste("Histogram of",names(physio[i])), xlab=names(physio[i]), las=1)
lines(density(physio[,i], bw=0.2))
rug(physio[,i])
}

## plot the empirical cumulative distribution
par(mfrow=c(3,4))
for(i in 1:ncol(neutral)){
plot(ecdf(neutral[,i]), do.points=FALSE, verticals=TRUE, main=paste("Cum Dist of",names(neutral[i])), xlab=names(neutral[i]), las=1)
}

par(mfrow=c(1,3))
for(i in 1:ncol(physio)){
plot(ecdf(physio[,i]), do.points=FALSE, verticals=TRUE, main=paste("Cum Dist of",names(physio[i])), xlab=names(physio[i]), las=1)
}

## Quantile-quantile (q-q)
par(mfrow=c(3,4))
for(i in 1:ncol(neutral)){
qqnorm(neutral[,i], main=paste("Q-Q plot of",names(neutral[i])), xlab=names(neutral[i]), las=1);qqline(neutral[,i])
}

par(mfrow=c(1,3))
for(i in 1:ncol(physio)){
qqnorm(physio[,i], main=paste("Q-Q plot of",names(physio[i])), xlab=names(physio[i]), las=1);qqline(physio[,i])
}


## Run some simulated data
x <- rt(250, df = 5)
qqnorm(x); qqline(x)
qqplot(qt(ppoints(250), df = 5), x, xlab = "Q-Q plot for t dsn")
qqline(x)

## Correlations (pearsons | spearman | kendall | chi square | Kruskal)
## Bivariate Relationship; 1 response 1 predictor
cat("\n\nHow much the variables are associated? ----- \nPositvie (evolve in a linear same direction) or \nnegative (evolve in a linear opposite direction) or \n0 (no linear association, but that doesnt mean there is no relation at all)", "\n")

## Association/Correlation Between Paired Samples (neutral | physio vs Developmental Stages)
b <- as.numeric(gl(4,3,12, label=c("0", "10", "21", "32")))			# development stages
Neutral_Spearman <- apply(neutral, 2, function(x) {cor.test(x, b, method="spearman")})
Physio_Spearman <- apply(neutral, 2, function(x) {cor.test(x, b, method="spearman")})

## Association/Correlation Between Paired Samples (physio vs neutral)
Spearman_neut_phy <- list()		## PI OR UI
z <- physio
for(i in 1:length(physio)){
Spearman_neut_phy[[i]] <- apply(neutral, 2, function(x) {cor.test(x, z[,i], method="spearman")})
}

## Difference in means among more than two groups
cat("\n\nAre the variables independent or not ?\n")
b <- gl(4,3,12, label=c("a.E", "b.V", "c.P", "d.J"))
Neutral_Kruskal <- apply(neutral, 2, function(x) {kruskal.test(x ~ b)})
Physio_Kruskal <- apply(physio, 2, function(x) {kruskal.test(x ~ b)})

## 1 predictor (2 levels), 1 dependent
## t-test
t.test(y, x)

### Test of Independence (non related variables)
## Chi squared
# Construct the contingency table; 1 respnse 1 predictor are nominal
tbl <- table(y,x)
chisq.test(tbl)
## 2-sample Wilcoxon (or Mann-Whitney)
## source : http://tinyurl.com/m2wb2b4
### END (diagnostics)


========================================
## Linear regression of neutral lipids
========================================

## GLM
## both glms can't work with logarithmic numbers nor NAs
# Fatty acids dataset (no trocophore)

## load datasets
setwd("C:/Dropbox/Workshop2013/Work/Paramètres Ecophysio/")
load("biometrics.Rdata", .GlobalEnv)
lsos(pat="")

## Preprocess dataset; Calculate mean FA content across stages of development 'egg inlcuded'
neutral
neutral2 <- fattyacids[-c(1:6),-c(2,16:17)]
neutral2
library(plyr)
means.neut <- ddply(neutral2, .(Treatment,Stage), numcolwise(mean))
means.neut

## LOAD GRowth and survival datasets for both treatments
## Growth = X axis
## Survival = Y axis
setwd("C:/Dropbox/Workshop2013/Work/R/Biometrie/")
biometrie <- read.table("Biometrie.txt", sep="\t", header=T)
names(biometrie)

## DATASET 3
## source : C:\Dropbox\Workshop2013\Work\Paramètres Ecophysio\Biometrics.xlsx
biometrics <- read.table("clipboard",sep="\t", header=T)
summary(biometrics)
names(biometrics)



## Datasets
sgrowth <- biometrie[, c(4,6,8)]
names(sgrowth)
survival <- biometrie[, c(4,6,9)]
names(survival)

summary(sgrowth)
summary(survival)
## Calculate mean Growth
means.growth <- ddply(sgrowth, .(Treatment,DPF.1), numcolwise(mean))
means.growth
gm.com <- rbind(c(0,50),				# size of the egg 50µm (?)a
			apply(subset(gm[gm$DPF.1 %in% c(13,15),]), 2, FUN=mean),
			apply(subset(gm[gm$DPF.1 %in% c(18,20),]), 2, FUN=mean),
			gm[gm$DPF.1==25,]); rownames(gm.com) <- stage

## Calculate mean Survival
survival <- read.table("survival.txt", sep="\t", header=T)
sm <- ddply(survival, .(DPF.1), numcolwise(mean))
sm.com <- rbind(sm[1,],				# size of the egg 5µm (?)
			apply(subset(sm[sm$DPF.1 %in% c(13,15),]), 2, FUN=mean),
			apply(subset(sm[sm$DPF.1 %in% c(18,20),]), 2, FUN=mean),
			sm[sm$DPF.1==25,]); rownames(sm.com) <- stage


####old

setwd("C:/Dropbox/Workshop2013/Work/Paramètres Ecopphysio/")
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

## Means of samples for the same stage
library(plyr)
samples <- as.numeric(gl(4,3,12, label=seq(1:4)))
neut_dd  <- data.frame(samples=samples, dataset)
dataset <- ddply(neut_dd, .(samples), numcolwise(mean))

stage <- c("Egg", "Veliger", "Pediveliger", "Juvenile")
rownames(dataset) <- stage ; dataset <- dataset[,-1]

## GRowth
library(plyr)
setwd("C:/Dropbox/Workshop2013/Work/R/datasets/")
stage <- c("Egg", "Veliger", "Pediveliger", "Juvenile")
growth <- read.table("Growth.txt", sep="\t", header=T)
gm <- ddply(growth, .(DPF.1), numcolwise(mean))
gm.com <- rbind(c(0,50),				# size of the egg 50µm (?)a
			apply(subset(gm[gm$DPF.1 %in% c(13,15),]), 2, FUN=mean),
			apply(subset(gm[gm$DPF.1 %in% c(18,20),]), 2, FUN=mean),
			gm[gm$DPF.1==25,]); rownames(gm.com) <- stage

## Survival
survival <- read.table("survival.txt", sep="\t", header=T)
sm <- ddply(survival, .(DPF.1), numcolwise(mean))
sm.com <- rbind(sm[1,],				# size of the egg 5µm (?)
			apply(subset(sm[sm$DPF.1 %in% c(13,15),]), 2, FUN=mean),
			apply(subset(sm[sm$DPF.1 %in% c(18,20),]), 2, FUN=mean),
			sm[sm$DPF.1==25,]); rownames(sm.com) <- stage

## Standardize x and y
x <- gm.com[,2]; y <- sm.com[,2]
ys <- (y-mean(y))/sd(y)
xs <- (x-mean(x))/sd(x)

# Transformation (1) double Square root
# Double square root regressions
srt.dataset <- dataset^0.25
d.sqrt <- NULL;var.expl <-NULL
for(j in 1:ncol(dataset)) {d.sqrt <- rbind(d.sqrt, lm(srt.dataset[,j]~ys+xs)$coefficients)
				var.expl <- sum(var.expl+var(lm(srt.dataset[,j]~ys+xs)$fitted.values))
				print(var(lm(srt.dataset[,j]~ys+xs)$fitted.values) / var(srt.dataset[,j]))
				}

cat("\n",var.expl / sum(apply(srt.dataset,2,var)),"\n")
summary(lm(srt.dataset~ys + xs))   ## This is where regressed coefficients can be found

plot(xs, ys, xlab = "x*(growth)", ylab = "y*(survival)", type = "n", asp = 1, cex.axis = 0.7, las=1)
text(xs, ys, labels = rownames(dataset), col = "forestgreen")
points(xs, ys, labels = rownames(dataset), col = "forestgreen", pch=17)
text(d.sqrt[,3:2], labels = colnames(dataset[,1:10]), col = "chocolate", font = 4, cex=0.7)
for(j in 1:ncol(dataset)) arrows(0, 0, 0.95*d.sqrt[j,3], 0.95*d.sqrt[j,2], col = "chocolate", angle = 15, length = 0.1)

