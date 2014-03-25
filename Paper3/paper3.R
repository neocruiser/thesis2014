##############################
#    Paper 3
##############################
source('c:/Dropbox/Workshop2013/Work/R/01funcs.R')
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\datasets\\weights\\")
#means <- read.table("meansAllPredictors.txt", sep="\t", header=T, row.names=3) ## 23 samples 50 K
#means <- read.table("means23samples15K.txt", sep="\t", header=T, row.names=3) ## 23 samples (all of them) 15K
means <- read.table("means22samples15K.txt", sep="\t", header=T, row.names=3) ## 22 samples 15K (without Juvenile Tiso replicate 1)
x <- t(means[,-c(1,2)])
dim(x)
## prepare features
y <- c(rep("E",3), rep("T",3), rep("VC",3),rep("PC",3),rep("JC",3),rep("VT",3),rep("PT",3), rep("JT",1))
y <- c(rep("E",3), rep("L",9), rep("PL",3), rep("L",6),rep("PL",1))
y <- c(rep("E",3), rep("L",9), rep("PL",3), rep("Lc",6),rep("PLc",1))
y <- c(rep("Healthy", 15), rep("Deficient",7))
y <- c(rep("None", 6), rep("Cocktail", 9), rep("Tiso",7))
y <- c(rep("E",3), rep("Lc",6), rep("PLc",6), rep("Lt",3),rep("PLt",4))
y <- c(rep("L",9), rep("PL",6), rep("L",3),rep("PL",4))
## prepare output (response varaible for only 22 samples JT1 has been excluded from analysis due to unadequate detection rules)
## prepare response variable (paper3)

##############################
# Filter subset selection
##############################
hi.x <- set.var(x,1,5000)
range(apply(hi.x, 2, var))
cor.x <- cor(hi.x); summary(cor.x[upper.tri(cor.x)])
## Unsupervised gene selection based on high variance (paper3)

require(mRMRe)
feature.select <- new("mRMRe.Data",data=data.frame(hi.x[,1:5000, drop=F])); gc()
## extract data
set.seed(1445612321)
locus.select <- new("mRMRe.Filter", data=feature.select, target_indices=1:5, levels=c(1000,1,1,1,1), continuous_estimator="spearman"); gc()
## feature select min redundant max relevant genes
range_correlation(locus.select,n=1,t=1,hi.x,method="spearman")	# n=gene(set w target_indices); t=different gene mashups
## view range of correlated selected features
locus <- locusRMR(locus.select,feature.select,n=1,t=1)
length(locus)
## extract LOCUS names
## Parallelized mRMR ensemble Feature selection

##############################
# Training/Testing
##############################
foo <- hi.x[,locus];dim(foo)
train <- sample(1:nrow(foo), 2*nrow(foo)/2.5); length(train)
test <- -train
y[test]; length(y[test])
## prepare training/testing sets
set.seed(1445612321)
dat <- data.frame(y=y, foo)
y <- as.vector(model.matrix(~y,dat)[,2])
#y <- c(rep(0,6),rep(1,9),rep(-1,7))	## dummy variables for 3 levels None, Coc, Tiso

##############################
# Feature extraction
##############################

## LASSO
grid <- 10^seq(10, -2, length=100)
require(glmnet)
set.seed(1445612321)
lasso.mod <- glmnet(foo[train,],y[train], alpha=1, lambda=grid, family = "gaussian", standardize=T, type.gaussian = "naive")
## (1) Model Training. Select the kernel function and associated kernel parameters
#coef(lasso.mod, s=0)
#par(mfrow = c(2,2));
plot(lasso.mod, xvar="lambda", label=T)
plot(lasso.mod, xvar="dev", label=T) 	## fraction deviance explained =R2
system.time(cv.out <- cv.glmnet(foo[train,], y[train], alpha=1, family="gaussian", standardize=T, nfolds = 10, type.gaussian="naive"))
## (2) Regularization. Cross validation for hyperparameter tuning. Optimization of model selection to avoid overfitting
par(mfrow = c(1,1)); plot(cv.out)
bestlam <- cv.out$lambda.min; bestlam
lasso.bestlam <- coef(lasso.mod, s=bestlam)
#bestlam <- cv.out$lambda.1se; bestlam
## Select the best hyperparameter
lasso.pred <- predict(lasso.mod, s=bestlam, newx=foo[test,], type="nonzero"); str(lasso.pred)
lasso.pred <- predict(lasso.mod, s=bestlam, newx=foo[test,], type="response"); lasso.pred
lasso.pred <- predict(lasso.mod, s=bestlam, newx=foo[test,], type="class")
table(lasso.pred, y[test])		## Confusion matrix for classification
mean((lasso.pred - y[test])^2)		## Test set MSE for regression
lasso.coef <- predict(lasso.mod, s=bestlam, type = "coefficients")
str(lasso.coef)
gene1leaf.mRMR <- lasso.coef	## SAVE to .Rdata
## show results
ind.lasso <- lasso.coef@i
loc.lasso <- lasso.coef@Dimnames[[1]]
final.select <- loc.lasso[ind.lasso]
foo <- row.names(means[rownames(means) %in% final.select,])
means[rownames(means) %in% final.select,2]
lasso.select <- x[,foo]
dim(lasso.select)
## extract genes from model selection
## LASSO (from the GLM package)

##############################
# Building classifier on subset model
##############################
## Choosing the right Hyper-parameters. GRID ANALYSIS
require(caret)
dat <- data.frame(y=y, lasso.select); dim(dat)
ctl=expand.grid(.size=seq(1,20,length=50), .decay=10^seq(-1,-5,length=50))	## nnet
#ctl=expand.grid(.C=, .sigma=0.27)	## svmRadial
#ctl=expand.grid(.C=10^seq(1,-4,length=100), .degree=10^seq(1,-1,length=40), .scale=10^seq(1,-3,length=40))	## svmPoly
#ctl=expand.grid(.alpha=seq(0.1,1,0.1), .lambda=10^seq(10, -4, length=100))	## glmnet
#ctl=expand.grid(.n.trees=100, .interaction.depth=1, .shrinkage=.001)	## Boosting (not working)
#ctl=expand.grid(.mtry=23)	## Random forest
set.seed(1445612321)
modelTune.reg(dat,train,test,method="nnet",folds=10,r=5,tune=10,ctl)	# Tune hyper-parameters
model.reg(dat,train,test,method="nnet",folds=10,r=5,tune=10)
## Regression
model.clas(dat,train,test,method="pls",folds=10,r=5,tune=10)
modelTune.clas(dat,train,test,method="pls",folds=10,r=5,tune=10,ctl)
## Classification
## Hyper-parameters tuning and model optimization

modelsRMSE <- read.table("clipboard", sep="\t", header=T);modelsRMSE
library(lattice)
xyplot(Gene1 + Gene2 + Gene3 + Gene4 + Gene5~Model, data=modelsRMSE, type=c("a","p"), pch =20,cex = 1,auto.key = list(space="top",points=T,lines=F))
xyplot(Tree1 + Tree2 + Tree3 + Tree4 + Tree5~Model, data=modelsRMSE, type=c("a","p"), pch =20,cex = 1,auto.key = list(space="top",points=T,lines=F))
xyplot(Gene1 + Gene2 + Gene3 + Gene4 + Gene5 ~ Model, data=modelsRMSE, type=c("a"),auto.key = list(space="top",points=F,lines=T))
## plot RMSE versus different learners

##############################
# Ensemble Learning (bagging)
##############################
require(caret)
require(foreach)
require(doSNOW)
cl <- makeCluster(3)
registerDoSNOW(cl)
#ctl=expand.grid(.size=0.5, .decay=0.03)	## nnet
ctl=expand.grid(.size=3, .decay=0)	## nnet
#ctl=expand.grid(.C=, .sigma=0.27)	## svmRadial
#ctl=expand.grid(.C=10^seq(1,-4,length=100), .degree=10^seq(1,-1,length=40), .scale=10^seq(1,-3,length=40))	## svmPoly
#ctl=expand.grid(.alpha=seq(0.1,1,0.1), .lambda=10^seq(10, -4, length=100))	## glmnet
#ctl=expand.grid(.n.trees=100, .interaction.depth=1, .shrinkage=.001)	## Boosting (not working)
#ctl=expand.grid(.mtry=23)	## Random forest
dat <- data.frame(y=y, lasso.select); dim(dat)
set.seed(1445612321)
bagging(dat[train,],dat[test,],m=1.1,ite=100,methods="nnet",tune=10)
## for testing
baggingTune(dat[train,],dat[test,],m=1.1,ite=100,methods="nnet",tune=10,gridZ=ctl)
## For tuning the hyper-parameters
bagging.clas(dat[train,],dat[test,],m=1.1,ite=100,methods="nnet",tune=10)
## For Classification
## END RUN
mean((modelbag - y[test])^2)		## Test set MSE for regression
stopCluster(cl)		## close cluster only after finishing w all modelse

## COMPUTE RMSE
ensemble.pred <- (gbbag+rfbag)/2
ensemble.pred <- (gbbag*2+rfbag)/3
ensemble.pred <- (gbbag+rfbag*2)/3
mean((ensemble.pred - y[test])^2)		## Test set MSE for regression
## Bagging

##############################
# Unsupervised PCR and Clustering
##############################
## START PCR
pr.out <- prcomp(lasso.select, scale=T)
par(mfrow = c(1,2))
y <- c(rep("E",3), rep("T",3), rep("VC",3),rep("PC",3),rep("JC",3),rep("VT",3),rep("PT",3), rep("JT",1))
plot(pr.out$x[,1:2],col=Cols(y),pch=19)
plot(pr.out$x[,c(1,3)],col=Cols(y),pch=19)
## 2D scatterplots
require(scatterplot3d)
par(mfrow = c(1,1))
scatterplot3d(pr.out$x[,1:3], pch=16, color=Cols(y), type="h", lty.hplot = "dashed", angle = 55, scale.y = .7)
## Plot 1 in 3D
require(rgl)
plot3d(pr.out$x[,1:3], size=5, col=Cols(y), type = "p", box=T, axes=T)
par3d(zoom = 1.1) # larger values make the image smaller
par3d(windowRect = c(50, 50, 500, 500)) # make the window large
text3d(pr.out$x[,1:3], text=y, font=1, cex=.6, adj=c(1,1))
## plot 2 in 3D
## START recording
M <- par3d("userMatrix") # save your settings to pass to the movie
movie3d(par3dinterp(userMatrix=list(M,rotate3d(M, pi, 1, 0, 0), rotate3d(M, pi, 0, 1, 0) ) ), duration = 5, fps = 50, movie = "MyMovie", dir = ("~/Downloads/PCA3D"))
## extract multiple frames for GIF implementation
#movie3d(spin3d(axis = c(0, 0, 1), rpm = 4), duration=15, movie="TestMovie", type="gif", dir = ("~/Downloads/PCA3D"))
#rgl.snapshot('PCA3D.png', fmt="png")
## extract 2
## END
## 3D scatterplot
summary(pr.out)
pve <- 100*pr.out$sdev^2/sum(pr.out$sdev^2)
plot(pr.out)
par(mfrow=c(1,2))
plot(pve, type="o", ylab = "PVE",xlab = "Principal Component", col="blue")
plot(cumsum(pve), type="o",ylab = "Cumultive PVE", xlab = "Principal Component", col="brown3")
## (Unsupervised) Principal Componenent analysis (reduced dimension) (paper3)

##############################
# Tryouts
##############################

cv.err <- rep(1:5)
for(i in 1:5){
glm.fit <- glm(mpg~poly(horsepower, i), data = Auto)
cv.err[i] <- cv.glm(Auto, glm.fit)$delta[1]
}
cv.err



pls_rf <- list(label = "PLS-RF",
               library = c("pls", "randomForest"),
               type = "Regression",
               ## Tune over both parameters at the same time
               parameters = data.frame(parameter = c('ncomp', 'mtry'),
                                       class = c("numeric", 'numeric'),
                                       label = c('#Components',
                                                 '#Randomly Selected Predictors')),
               grid = function(x, y, len = NULL) {
                 grid <- expand.grid(ncomp = seq(1, min(ncol(x) - 1, len), by = 1),
                                     mtry = 1:len)
                 ## We can't have mtry > ncomp
                 grid <- subset(grid, mtry <= ncomp)
               },
               loop = NULL,
               fit = function(x, y, wts, param, lev, last, classProbs, ...) {
                 ## First fit the pls model, generate the training set scores,
                 ## then attach what is needed to the random forest object to
                 ## be used later
                 ## plsr only has a formula interface so create one data frame
                 dat <- x
                 dat$y <- y
                 pre <- plsr(y~ ., data = dat, ncomp = param$ncomp)
                 scores <- predict(pre, x, type = "scores")
                 colnames(scores) <- paste("score", 1:param$ncomp, sep = "")
                 mod <- randomForest(scores, y, mtry = param$mtry, ...)
                 mod$projection <- pre$projection
                 mod
               },
               predict = function(modelFit, newdata, submodels = NULL) {
                 ## Now apply the same scaling to the new samples
                 scores <- as.matrix(newdata)  %*% modelFit$projection
                 colnames(scores) <- paste("score", 1:ncol(scores), sep = "")
                 ## Predict the random forest model
                 predict(modelFit, scores)
               },
               prob = NULL,
               varImp = NULL,
               predictors = function(x, ...) rownames(x$projection),
               levels = function(x) x$obsLevels,
               sort = function(x) x[order(x[,1]),])

train <- sample(1:nrow(foo), 2*nrow(foo)/2.4); length(train)
test <- -train
y[test]; length(y[test])
y <- c(rep("L",9), rep("PL",6), rep("L",3),rep("PL",4))
dat <- data.frame(y=y, lasso.select); dim(dat)
require(nnet)
stages <- class.ind(dat[,1])
dat <- data.frame(y=stages[,1], lasso.select); dim(dat)
## create continuous regression matrix
metaCtrl <- trainControl(method = "repeatedcv", repeats=3)
set.seed(10934)
Profile <- train(dat[,-1], dat[,1], method=pls_rf, preProc=c("center","scale"),tuneLength=10,ntree=500,trControl=metaCtrl)
## Customized Train function with regression model fitting (Not working)

require(gbm)
set.seed(135662)
boosted.model <- gbm(y~., data=dat[train,], distribution = "gaussian",n.trees=1000,interaction.depth = 4)
## BOOSTING (sample size too small)
require(randomForest)
bagging.model <- randomForest(y~., data=dat, subset=train, mtry=23,importance=T)
bagging.model
bagging.pred <- predict(bagging.model, newdata = dat[test,])
mean((bagging.pred-dat[test,1])^2)
plot(bagging.pred,dat[test,1],type="b")
## BAGGING (not working)


sd.data <- scale(lasso.select)
par(mfrow = c(1,1))
data.dist <- dist(sd.data)
plot(hclust(data.dist), labels=y, main="Complete Linkage")
plot(hclust(data.dist, method="single"), labels=y, main="Complete Linkage")
plot(hclust(data.dist, method="average"), labels=y, main="Average Linkage")
abline(h=18, col="red")
## Hierarchicla clustering
hc.out <- hclust(data.dist)
hc.clusters <- cutree(hc.out, 4)
table(hc.clusters, y)
hc.out
## cut the dendrogram at the height that yield a particular number
set.seed(1)
km.out <- kmeans(sd.data, 4, nstart = 20)
km.clusters <- km.out$cluster
table(km.clusters, hc.clusters)
## K means clustering
hc.out <- hclust(dist(pr.out$x[,1:6]))
plot(hc.out, labels=y, main="Hier. Clust., on First Six Score Vectors")
table(cutree(hc.out, 4), y)
## perform clustering on the first 6 rpincipal components

set.seed(2)
train <- sample(1:nrow(x), 2*nrow(x)/2.5)
test <- -train
dat <- data.frame(x[train,1:5000], y=as.factor(y[train]))
require(e1071)
out <- svm(y~., data=dat, kernel = "radial",cost=10, gamma=1)
summary(out)
table(out$fitted, dat$y)
## Training observations
dat.te <- data.frame(x[test,1:5000], y=as.factor(y[test]))
pred.te <- predict(out, newdata=dat.te)
table(pred.te, dat.te$y)
cat("\nMisclassification Error rate",(4/22)*100,"\n")
## support vector regression

require(leaps)
train <- sample(1:43000, 1000)
dat <- data.frame(x[,train], y=as.factor(y))
regfit.fwd <- regsubsets(y~., data=dat, method="backward", nvmax = 1000)
reg.summary <- summary(regfit.fwd)
plot(reg.summary$rss, type="l")
plot(reg.summary$adjr2, type="l")
## backward selection (not working)


train <- sample(1:nrow(foo), 2*nrow(foo)/2.5); length(train)
test <- -train
y[test]; length(y[test])
y <- c(rep("L",9), rep("PL",6), rep("L",3),rep("PL",4))
dat <- data.frame(y=y, lasso.select); dim(dat)
require(nnet)
stages <- class.ind(dat[,1])
dat <- data.frame(y=stages[,1], lasso.select); dim(dat)
## create continuous regression matrix
training <- dat[train,]
testing <- dat[-train,]
ctl=expand.grid(.mtry=17)	## Random forest
ctl=expand.grid(.fL=1, .usekernel=1)	## Random forest
ctl=expand.grid(.alpha=0.1, .lambda=0.9)	## Glmnet (lasso)
set.seed(156934)
## START RUN
system.time(
glmnetbag <- bagging(training,testing,m=3,ite=200,meth="glmnet",gridZ=ctl,tune=5))
## END RUN
sqrt((sum((stages[test,1] - Predd)^2))/nrow(stages[test,]))		## compute RMSE (REGRESSION)
## BAGGING CUSTOM FUNCTION (not working)



##############################
# SAVE
##############################

setwd("C:/Dropbox/Workshop2013/Work/R/ANN")
lsos(pat="locus.select|*.mRMR")
save(list=ls(pattern="*.mRMR"),file="lassoSelected.Rdata")	## save
save(list=ls(pattern="locus.select"),file="mRMRselected.Rdata")	## save
load("lassoSelected.Rdata", .GlobalEnv)
