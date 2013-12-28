## IMPORTANT***
## Load the summarySE function from 01funcs

=========================================================================
#	SECTION --								#
=========================================================================
========================
  Physiological studies
========================
## in the case of repeated measurements on the same group (e.g., of animals, individuals, cultures, or reactions),
## CIs or SE bars are irrelevant to comparisons within the same group
## Read biometrie file
setwd("C:/Dropbox/Workshop2013/Work/R/Biometrie/")
biometrie <- read.table("Biometrie.txt", sep="\t", header=T)
names(biometrie)
library(ggplot2)

  ## GROWTH
growth <- biometrie[1:695,c(6,8)]
growth[1:20,]
growth <- growth[complete.cases(growth),]
		#write.table(growth, "growth.txt", sep="\t")
dtfm <- summarySE(growth, measurevar="Growth", groupvars="DPF.1")
# Confidence interval of the mean of GROWTH
ggplot(dtfm, aes(x=DPF.1, y=Growth)) +
  geom_errorbar(aes(ymin=Growth-se, ymax=Growth+se), width=.4) +
  geom_line() +
  geom_point(size=3.5, shape=21, fill="white")


  ## MORTALITY
names(biometrie)
survival <- biometrie[1:694, c(6,9)]
names(survival)
?complete.cases
survival <- survival[complete.cases(survival),]   # Remove NA rows
		# write.table(survival, "survival.txt", sep="\t")
s.dft <- summarySE(survival, measurevar="Mortality", groupvars="DPF.1")
# Confidence interval of the mean of Mortality
ggplot(s.df, aes(x=DPF.1, y=Mortality)) +
  geom_errorbar(aes(ymin=Mortality-se, ymax=Mortality+se), width=.4) +
  geom_line() +
  geom_point(size=3.5, shape=22, fill="white")


  ## PERFORMANCE
performance <- biometrie[1:694, c(6,10)]
performance <- performance[complete.cases(performance),]
performance
		#write.table(performance, "performance.txt", sep="\t")
p.df <- summarySE(performance, measurevar="Performance", groupvars="DPF.1")
# Confidence interval of the mean of Performance
ggplot(p.df, aes(x=DPF.1, y=Performance)) +
  geom_errorbar(aes(ymin=Performance-sd, ymax=Performance+sd), width=.4) +
  geom_line() +
  geom_point(size=3.5, shape=21, fill="white")


  ## COMPETENCE
competence <- biometrie[1:694, c(6,11)]
competence <- competence[complete.cases(competence),]
competence
		#write.table(competence, "competence.txt", sep="\t")
c.df <- summarySE(competence, measurevar="Compétence", groupvars="DPF.1")
# Confidence interval of the mean of Competence
ggplot(c.df, aes(x=DPF.1, y=Compétence)) +
  geom_errorbar(aes(ymin=Compétence-sd, ymax=Compétence+sd), width=.4) +
  geom_line() +
  geom_point(size=3.5, shape=21, fill="white")



## Quality
quality <- biometrie[1:694, c(6,12)]
quality <- quality[complete.cases(quality),]
quality
		#write.table(quality, "quality.txt", sep="\t")
q.df <- summarySE(quality, measurevar="Condition", groupvars="DPF.1")
# Confidence interval of the mean of Condition
ggplot(q.df, aes(x=DPF.1, y=Condition)) +
  geom_errorbar(aes(ymin=Condition-sd, ymax=Condition+sd), width=.7) +
  geom_line() +
  geom_point(size=3.5, shape=21, fill="white")



=========================================================================
#	SECTION --								#
=========================================================================

## Anova of indices & lipids
## + Diagnostics
## Profil Lipides 2011.xlsx
## PHYSIO*********************
library(ggplot2)
samples <- gl(4,3,12, label=c("a.E", "b.V", "c.P", "d.J"))
summary(samples)
indices <- read.table("C:\\Dropbox\\Workshop2013\\Work\\Paramètres Ecophysio\\indices.txt", sep="\t", header=T, skip=1, row.names=1)
physio <- data.frame(PI=indices[-c(7:9,13:15,19:21),7],
					OR=indices[-c(7:9,13:15,19:21),11],
					UI=indices[-c(7:9,13:15,19:21),20])

## NEUTRAL**************
## Fatty acids dataset (********no trocophore, no sums*********)
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



## Diagnostics (Anova, Multiple comparison)
## one way ANOVA (2 or more conditions)
## START (run)
cat("\n R2 represents the % of variance in the response variable accounted for or explained by the variability in the predictor variable\n\n")

## Anova & Tukey HSD multicomparison test
a <- gl(4,3,12, label=c("a.E", "b.V", "c.P", "d.J"))

Neutral_anova <- apply(neutral, 2, function(x) {
list(
fit <- aov(x ~ a),
p <- drop1(fit, ~., test="F"),
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),
TukeyHSD(fit, ordered=T)
)})

Physio_anova <- apply(physio, 2, function(x) {
list(
fit <- aov(x ~ a),
p <- drop1(fit, ~., test="F"),
print(c("R^2=",p[2,2]/p[2,3])),
summary(fit),
TukeyHSD(fit, ordered=T)
)})

## Bonferroni correction (multiple comparison)
Neutral_nonMultipleComparison <- apply(neutral, 2, function(x) pairwise.t.test(x, a, p.adj= "none"))
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
## END (run)

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


## Plot
## Either physio or neutral datasets
samples <- as.numeric(gl(4,3,12, label=seq(1:4)))
dtfm <- summarySE(physio, measurevar="OR", groupvars="samples")

## colname of X (must be the same as ^^^^)
plot(dtfm$PI, type="n",las=1, ylim=c(0,140),xlab="Stages", ylab="Mol %", yaxt="n" , xaxt="n")
axis(side=1, at=c(1:4), labels=c("E","V","P","J"))
axis(2, at=c(seq(1,15,5), seq(60,170,20)), las=1)
axis(4, at=dtfm$OR,labels=round(dtfm$OR, digits=2),las=2, cex.axis=0.7, tck=-.01)
legend(x="top", c("PI", "UI", "OR"), pch=c(21,22,24), lty=1, bty="n", pt.bg=c("chocolate"))

## SET X for each indice
x <- dtfm$OR
lines(x, col="black")	; se <- dtfm$se
arrows(dtfm$samples,x, dtfm$samples,x+se, code=2, angle=90, length=0.1, col="black")
points(x, pch=24, type="p", bg="chocolate", cex=1.3)
	## DRAW a gap plot
	# require(plotrix)
	# gap.plot(dtfm$PI, gap=c(21,60), ytics=c(5,15,seq(80,170,20)), type="n",las=1, ylim=c(0,170), xtics=c(1:4), xticlab=c("E","V","P","J"))

=========================================================================
#	SECTION --								#
=========================================================================
## Linear regression of neutral lipids
========================================

## Work with 3 replicates COCktail
## Fatty acids dataset (********no trocophore, no sums*********)
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

## Boxplot for all neutral variables
layout(matrix(c(1:10),2,5))
for(i in 1:ncol(neutral)) boxplot.n(neutral[,i] ~ samples, shrink=0.8, notch=F, boxwex=0.25, las=1)

## Barplot for all neutral variables (means +sd)
library(plyr)
samples <- as.numeric(gl(4,3,12, label=seq(1:4)))
#logs.n <- log(neutral)
neut_dd  <- data.frame(samples=samples, neutral)
fit <- ddply(neut_dd, .(samples), numcolwise(mean))
barplot(as.matrix(t(fit[,-1])),
		col=gray.colors(10),
		beside=T, names.arg=c("Egg","Veliger","Pediveliger","Juvenile"),
		las=1, border=F,
		#legend.text=colnames(neutral),args.legend = list(x = "topleft")
		)

## Plot with base plotting
## Redefine samples as to DRAW connecting lines
samples <- as.numeric(gl(4,3,12, label=seq(1:4)))
dtfm <- summarySE(neutral, measurevar="TAG", groupvars="samples")									# Change variable (3)
## colname of X (must be the same as ^^^^)
plot(dtfm$TAG, type="n",las=1, ylim=c(0,60), 														# change ylim (1)
xlab="Stages", ylab="Mol %", , xaxt="n")
abline(h=c(31,4.2), col=c("chocolate", "forestgreen"))												# change "h" dietary coordinates (2)
axis(side=1, at=c(1:4), labels=c("E","V","P","J"))
legend(x="top", c("PL", "TAG", "AMPL","ST"), pch=c(21,22,23,24), lty=1, bty="n",
pt.bg=c("chocolate","chocolate","chocolate","forestgreen")) 										# change legend (4)

## loop for every fraction
x <- dtfm$TAG																						# Change variable (5)
lines(x, col="black")	; se <- dtfm$se
arrows(dtfm$samples,x, dtfm$samples,x+se, code=2, angle=90, length=0.1, col="black")
points(x, pch=22, type="p", bg="forestgreen", cex=1.3)												# Change pch (6) & color (7)

# arrows(dtfm$samples,x-se, dtfm$samples,x+se, code=3, angle=90, length=0.1, col="chocolate")a


#	Other example for drawing error bars (package Hmisc)
#	errbar(dtfm$samples, dtfm$DHA, dtfm$DHA+se, dtfm$DHA-se, add=T)

#########################################################################
#	SECTION -4-								#
#########################################################################
## GLM
## both glms can't work with logarithmic numbers nor NAs
# Fatty acids dataset (no trocophore)

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

