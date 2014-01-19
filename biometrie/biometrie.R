## IMPORTANT***
## Load the summarySE function from 01funcs

========================
  Physiological studies
========================

## load datasets
setwd("C:/Dropbox/Workshop2013/Work/Paramètres Ecophysio/")
load("biometrics.Rdata", .GlobalEnv)
lsos(pat="")
windows(record=T)
library(ggplot2)

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
names(fattyacids)
neutral <- fattyacids[-c(1:9),-c(1:2,5,20:21)]
neutral

## DATASET 3
## source : C:\Dropbox\Workshop2013\Work\Paramètres Ecophysio\Biometrics.xlsx
biometrics <- read.table("clipboard",sep="\t", header=T)
summary(biometrics)
names(biometrics)

## GROWTH for both diets COC and TISO
growth <- biometrics[biometrics$Variable == "Growth",]
names(growth)
growth <- growth[complete.cases(growth),]   # Remove NA rows
DF <- summarySE(growth[-c(1:3),], measurevar="Physio", groupvars=c("DPF","Treatment"))
DF
## order the x axis ascending (very important)
DF$DPF <- factor(DF$DPF, levels=c("6","8","11","13","15","18","20","22","25"))
# CONFIDENCE INTERVAL OF THE MEAN OF GROWTH
ggplot(DF, aes(x=DPF, y=Physio)) +
  geom_errorbar(aes(ymin=Physio-se, ymax=Physio+se), width=.4) +
  geom_line(aes(group=Treatment)) +
  geom_point(size=4.4, fill="black",aes(shape=Treatment)) +
    theme_bw()

## Mortality COC and TISO
survival <- biometrics[biometrics$Variable == "Mortality",]   # Remove NA rows
head(survival)
survival <- survival[complete.cases(survival),]   # Remove NA rows
DF <- summarySE(survival, measurevar="Physio", groupvars=c("DPF","Treatment"))
DF$DPF <- factor(DF$DPF, levels=c("1","4","6","8","11","13","15","18","20","22","25"))
# Confidence interval of the mean of Mortality
ggplot(DF, aes(x=DPF, y=Physio)) +
  geom_errorbar(aes(ymin=Physio-se, ymax=Physio+se), width=.4) +
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
lsos(pat=".*indices|fatty.*|neutral|biometrics")
save(list=ls(pattern="indices|fatty.*|neutral|biometrics"),file="biometrics.Rdata")


====================
    Factorial MANOVA
====================

## source : http://tinyurl.com/y93l9xq
    ## 2x2 factorial MANOVA
    ## load datasets
setwd("C:/Dropbox/Workshop2013/Work/Paramètres Ecophysio/")
load("biometrics.Rdata", .GlobalEnv)
lsos(pat=".*indices|fatty.*|neutral|biometrics")

Y <- as.matrix(neutral[,-c(1:2)])
## transform neutral to a numeric matrix

Diet <- gl(2,9,18, labels=c("COC","Tiso"))
Time <- gl(3,3,18, labels=c("Veliger","Pediveliger","Juvenile"))
## create factors for Treatment and Stage of development

X <- as.matrix(neutral[,-c(1:2)])
neutral.trans <- asin(sqrt(X/100))
neutral.trans
## Arcsine transformation of the data to achive constant variance
Y <- as.matrix(neutral.trans[,c(6:7,16)])
fit <- manova(Y ~ Diet*Time)				## 2x2 factorial MANOVA
options(digits=3)
summary.manova(fit, test="W")				## summary on the pairs of coordinates
summary.aov(fit)				##summary on the individual coordinates (type I SS)
## generate summary for ANOVA single entries or MANOVA overall entries (change neutral.trans colmns)

neutral.tukey <- apply(Y, 2, function(x) {
list(
fit <- aov(x ~ diet*samples),
p <- drop1(fit, ~., test="F"),		## type III SS and F test
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),				## type I SS Anova table
TukeyHSD(fit, ordered=T)
)})
## Multiple comparisons


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

lsos(pat=".*indices|fatty.*|neutral|biometrics")
## show DATASETS (created from clipboard, found in the Paramètres Ecophysio (/work directory))

## All samples + Both treatments
anova.neutral.all <- apply(neutral[,-c(1:2)], 2, function(x) {
list(
fit <- aov(x ~ diet*samples),
p <- drop1(fit, ~., test="F"),		## type III SS and F test
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),				## type I SS Anova table
TukeyHSD(fit, ordered=T)
)})

## All samples + one treatment
## COC treatment
anova.neutral.between <- apply(neutral[,-c(1:2)], 2, function(x) {
list(
fit <- aov(x ~ samples + diet),
p <- drop1(fit, ~., test="F"),
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),
TukeyHSD(fit, ordered=T)
)})

## Tiso treatment
anova.neutral.single.Tiso <- apply(neutral[neutral$Treatment=="Tiso",-c(1:2)], 2, function(x) {
list(
fit <- aov(x ~ samples[1:9]),
p <- drop1(fit, ~., test="F"),
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),
TukeyHSD(fit, ordered=T)
)})

## Cocktail treatment
anova.neutral.single.coc <- apply(neutral[neutral$Treatment=="Cocktail",-c(1:2)], 2, function(x) {
list(
fit <- aov(x ~ samples[1:9]),
p <- drop1(fit, ~., test="F"),
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),
TukeyHSD(fit, ordered=T)
)})


====================
    Additional tools
====================

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
lsos(pat=".*indices|fatty.*|neutral|biometrics")

## Preprocess dataset; Calculate mean FA content across stages of development 'egg inlcuded'
neutral
neutral2 <- fattyacids[-c(1:6),-c(2,16:17)]
neutral2
library(plyr)
means.neut <- ddply(neutral2, .(Treatment,Stage), numcolwise(mean))
means.neut

## Calculate mean Growth
growth <- biometrics[biometrics$Variable == "Growth",]
head(growth)
means.growth <- ddply(growth, .(Treatment,Stage), numcolwise(mean))
means.growth

## Calculate mean Survival
survival <- biometrics[biometrics$Variable == "Mortality",]
head(survival)
means.survival <- ddply(survival, .(Treatment, Stage), numcolwise(mean))
means.survival

## Standardize x and y
x <- means.growth[c(5,2:4),4]
y <- means.survival[c(5,2:4),4]
ys <- (y-mean(y))/sd(y)
xs <- (x-mean(x))/sd(x)

# Transformation (1) double Square root
# Double square root regressions
srt.dataset <- means.neut[c(4,1:3),-c(1:3)]^0.25
d.sqrt <- NULL
var.expl <-NULL
for(j in 1:(dim(means.neut)[2]-3)) {
    d.sqrt <- rbind(d.sqrt, lm(srt.dataset[,j]~ys+xs)$coefficients)
    var.expl <- sum(var.expl+var(lm(srt.dataset[,j]~ys+xs)$fitted.values))
    print(var(lm(srt.dataset[,j]~ys+xs)$fitted.values) / var(srt.dataset[,j]))
	}

cat("\n",var.expl / sum(apply(srt.dataset,2,var)),"\n")
summary(lm(srt.dataset~ys + xs))   ## This is where regressed coefficients can be found

plot(xs, ys, xlab = "x*(growth)", ylab = "y*(survival)", type = "n", asp = 1, cex.axis = 0.7, las=1)
means.neut$Stage <- factor(means.neut$Stage, levels=c("Egg","Veliger","Pediveliger","Juvenile"))
text(xs, ys, labels = means.neut[c(4,3:1),2], col = "forestgreen")
points(xs, ys, col = "forestgreen", pch=17, cex=1.5)
text(d.sqrt[,3:2], labels = colnames(means.neut[,-c(1:3)]), col = "chocolate", font = 4, cex=0.7)
for(j in 1:(dim(means.neut)[2]-3)) arrows(0, 0, 0.95*d.sqrt[j,3], 0.95*d.sqrt[j,2], col = "chocolate", angle = 15, length = 0.1)



## Tiso treatment
## Standardize x and y
x <- means.growth[c(5,6:8),4]
y <- means.survival[c(5,6:8),4]
ys <- (y-mean(y))/sd(y)
xs <- (x-mean(x))/sd(x)

# Transformation (1) double Square root
# Double square root regressions
srt.dataset <- means.neut[c(4,5:7),-c(1:3)]^0.25
d.sqrt <- NULL
var.expl <-NULL
for(j in 1:(dim(means.neut)[2]-3)) {
    d.sqrt <- rbind(d.sqrt, lm(srt.dataset[,j]~ys+xs)$coefficients)
    var.expl <- sum(var.expl+var(lm(srt.dataset[,j]~ys+xs)$fitted.values))
    print(var(lm(srt.dataset[,j]~ys+xs)$fitted.values) / var(srt.dataset[,j]))
	}

cat("\n",var.expl / sum(apply(srt.dataset,2,var)),"\n")
summary(lm(srt.dataset~ys + xs))   ## This is where regressed coefficients can be found

plot(xs, ys, xlab = "x*(growth)", ylab = "y*(survival)", type = "n", asp = 1, cex.axis = 0.7, las=1)
means.neut$Stage <- factor(means.neut$Stage, levels=c("Egg","Veliger","Pediveliger","Juvenile"))
text(xs, ys, labels = means.neut[c(4,7:5),2], col = "forestgreen")
points(xs, ys, col = "forestgreen", pch=17, cex=1.5)
text(d.sqrt[,3:2], labels = colnames(means.neut[,-c(1:3)]), col = "chocolate", font = 4, cex=0.7)
for(j in 1:(dim(means.neut )[2]-3)) arrows(0, 0, 0.95*d.sqrt[j,3], 0.95*d.sqrt[j,2], col = "chocolate", angle = 15, length = 0.1)


