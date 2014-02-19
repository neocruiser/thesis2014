
=========================================================================
#	SECTION --								#
=========================================================================
##### USING NEURALNET PACKAGE GENERATE A PLOT OF THE NEURAL NETWORK
require(ISLR)
data(Smarket)
head(Smarket)
dim(Smarket)

# source : http://tinyurl.com/oekteno
require(neuralnet)
mod1 <- neuralnet(Year~Lag2+Lag3+Lag4+Lag5+Volume+Today, data=Smarket, hidden=10)

#import the function from Github
library(devtools)
plot.nnet(mod1)
## downlaod and install function plot.nnet
## this function works with nnet, neuralnet, and RSNNS packages
source_url('https://gist.github.com/fawda123/7471137/raw/cd6e6a0b0bdb4e065c597e52165e5ac887f5fe95/nnet_plot_update.r')

=========================================================================
#	SECTION --								#
=========================================================================
## RECURSIVE FEATURE ELIMINATION

## DATA
library(caret)
library(nnet)
data(Caravan)
data(iris)
irisdata <- read.csv("http://www.heatonresearch.com/dload/data/iris.csv",head=TRUE,sep=",")
data(Smarket, package="ISLR") ## binary categorical


set.seed(34)
## classification only
nbFuncs$summary <- twoClassSummary
treebagFuncs$summary <- twoClassSummary
rfFuncs$summary <- twoClassSummary
gamFuncs$summary <- twoClassSummary
ctrl <- rfeControl(functions=treebagFuncs, method="repeatedcv", repeats=3, verbose=F)
Profile <- rfe(Caravan[,-c(1:30,50:60,86)], Caravan[,86], sizes=c(20,30,40,50), rfeControl=ctrl, metric="ROC")

## Regression only
## create dummy dependent variables
x=Smarket
NNi <- class.ind(x[,dim(x)[2]])
NN <- data.frame(x[,-dim(x)[2]], species=NNi[,1])
## compute backward selection
ctrl <- rfeControl(functions=caretFuncs, method="svmRadial", verbose=F, fit=FALSE)
Profile <- rfe(NN[,-dim(x)[2]], NN[,dim(x)[2]], sizes=c(2,3,4,5,8), rfeControl=ctrl)

Profile
xyplot(Profile$results$ROC ~ Profile$results$Variables, type=c("p","g","l"),auto.key=T)
Profile$fit
head(Profile$resample)
plot(Profile,type="b")


====================
    Visualization of the dataset
====================

library(AppliedPredictiveModeling)
transparentTheme(trans=.4)
## activate plotting theme

head(iris)
featurePlot(x = iris[,-5],
            y = iris$Species,
            plot = "ellipse",
            auto.key = list(columns=3))
## scatterplots (source : http://tinyurl.com/k6m6m5r)

=========================================================================
#								SECTION --								#
=========================================================================

## COMPARE DIFFERENT CLASSIFIERS AND TUNING
setwd("C:/Dropbox/Workshop2013/Work/R/ANN/")
load("caret.Rdata", .GlobalEnv)
lsos(pat="Train.*")

require(caret)
require(nnet)
require(kernlab)
require(caTools)

#### CROSS VALIDATION ON THE LEARNER
## DATA
data(Smarket, package="ISLR") ## binary categorical
data(Caravan, package="kernlab")
data(spam, package="kernlab")
NN=spam[1:2000,]

## create binary dataset
x = market				## original dataset
## choose this "for classification"
NN <- data.frame(x[,-dim(x)[2]], species=x[,dim(x)[2]])		## all dataset categorical (class)
# OR this "for regression"
NNi <- class.ind(x[,dim(x)[2]])		## create continuous variable for regression
NN <- data.frame(x[,-dim(x)[2]], species=NNi[,1])	## all dataset binary continuous
# finish
y = as.matrix(NN[,dim(NN)[2]]) 		## only species

## CREATE CROSS VALIDATION K FOLDS FOR TRAINING
trainList <- createDataPartition(y, p = 0.8, list = FALSE)
str(trainList)
training <- NN[trainList,]
testing <- NN[-trainList,]

#### SVM Radial (D)
set.seed(21)
#training
svmTrain <- train(species~.,
			data=training,
			method="svmRadial",
			trControl=trainControl(method="repeatedcv",	number=10, repeats=5
				#classProbs=T,
				#summaryFunction=twoClassSummary		## ROC curve
				),
			tuneGrid=data.frame(.C=2.2, .sigma=0.27),
			preProc=c("center","scale"),
			tuneLength=5
			#metric="ROC"
			)
svmTrain

plot(svmTrain)

## RMSE + R2
svmTest <- predict(svmTrain, newdata=testing)
cat('\nComputeRMSE','\nRMSE=',
	sqrt((sum((testing$species-svmTest)^2))/nrow(testing)),
	'\n',
	'\nCompute Rsquared', '\nR2=',
	svmTrain$results$Rsquared,
	'\n')

## For binary classification (only)
Predd <- predict(svmTrain, newdata=testing, type="prob")
colAUC(Predd, testing, plot=T)

## COMPUTE A CONFUSION MATRIX
confusionMatrix(data=svmTest, testing$species)




#### SVM Radial weight (C)
set.seed(34)
#training
srwTrain <- train(species~.,
			data=training,
			method="svmRadialWeights",
			trControl=trainControl(method="repeatedcv",	number=10, repeats=5
				#classProbs=T,
				#summaryFunction=twoClassSummary		## ROC curve
				),
			tuneGrid=data.frame(.C=4, .sigma=0.16,.Weight=1),
			preProc=c("center","scale"),
			#metric="ROC",
			tuneLength=5
			)
srwTrain

plot(srwTrain)
## RMSE + R2
srwTest <- predict(srwTrain, newdata=testing)
cat('\nComputeRMSE......','RMSE=',
	sqrt((sum((testing$species-srwTest)^2))/nrow(testing)),
	'\n',
	'\nCompute Rsquared.......', 'R2=',
	srwTrain$results$Rsquared
	)

## For binary classification (only)
## AUC
srwTest <- predict(srwTrain, newdata=testing)
cat('\nComputeROC........','ROC=',
	srwTrain$results$ROC,
	'\n')


Predd <- predict(srwTrain, newdata=testing, type="prob")
colAUC(Predd, testing, plot=T)

## COMPUTE A CONFUSION MATRIX
confusionMatrix(data=srwTest, testing$species)



#### nnet (D)
set.seed(21)
#training
nnetTrain <- train(species~.,
				data=training,
				method="nnet",
				tuneGrid=data.frame(.size=1, .decay=5e-6),
				trControl=trainControl(method="repeatedcv",	number=10, repeats=5
					#classProbs=T,
					#summaryFunction=twoClassSummary		## ROC curve
					),
				#metric="ROC",
				preProc=c("center","scale")
				)
nnetTrain

plot(nnetTrain)

## RMSE + R2
nnetTest <- predict(nnetTrain, newdata=testing)
cat('\nComputeRMSE','\nRMSE=',
	sqrt((sum((testing$species-nnetTest)^2))/nrow(testing)),
	'\n',
	'\nCompute Rsquared', '\nR2=',
	nnetTrain$results$Rsquared,
	'\n')

## For binary classification (only)
## AUC
nnetTest <- predict(nnetTrain, newdata=testing)
cat('\nComputeROC........','ROC=',
	nnetTrain$results$ROC,
	'\n')


Predd <- predict(nnetTrain, newdata=testing, type="prob")
colAUC(Predd, testing, plot=T)

## COMPUTE A CONFUSION MATRIX
confusionMatrix(data=nnetTest, testing$species)



#### Random forest ((D)) ---- Bagging is already an integrated feature
set.seed(21)
#training
rfTrain <- train(species~.,
			data=training,
			method="rf",
			tuneGrid=data.frame(.mtry=(17)),
			trControl=trainControl(method="repeatedcv",	number=10, repeats=5
					#classProbs=T,
					#summaryFunction=twoClassSummary		## ROC curve
					),
			preProc=c("center","scale")
			#metric="ROC"
			)
rfTrain

plot(rfTrain)

## RMSE + R2
rfTest <- predict(rfTrain, newdata=testing)
cat('\nComputeRMSE','\nRMSE=',
	sqrt((sum((testing$species-rfTest)^2))/nrow(testing)),
	'\n',
	'\nCompute Rsquared', '\nR2=',
	rfTrain$results$Rsquared,
	'\n')

#### For binary classification (only)
## AUC
rfTest <- predict(rfTrain, newdata=testing)
cat('\nComputeROC........','ROC=',
	rfTrain$results$ROC,
	'\n')

Predd <- predict(rfTrain, newdata=testing, type="prob")
colAUC(Predd, testing, plot=T)

## COMPUTE A CONFUSION MATRIX
confusionMatrix(data=rfTest, testing$species)





#### Boosted trees ((D)) ---- Generalized Boosted Regression Model
set.seed(21)
#training
gbmTrain <- train(species~.,
			data=training,
			method="gbm",
			tuneGrid=data.frame(.interaction.depth=c(1,5,9), .n.trees=(1:15)*100,.shrinkage=0.1),
			trControl=trainControl(method="repeatedcv",	number=10, repeats=5
					#classProbs=T,
					#summaryFunction=twoClassSummary		## ROC curve
					),
			preProc=c("center","scale")
			#metric="ROC"
			)
gbmTrain

plot(gbmTrain)
gbmTrain$finalModel

## RMSE + R2
gbmTest <- predict(gbmTrain, newdata=testing)
cat('\nComputeRMSE','\nRMSE=',
	sqrt((sum((testing$species-gbmTest)^2))/nrow(testing)),
	'\n',
	'\nCompute Rsquared', '\nR2=',
	gbmTrain$results$Rsquared,
	'\n')

### For binary classification (only)
## AUC
gbmTest <- predict(gbmTrain, newdata=testing)
cat('\nComputeROC........','ROC=',
	gbmTrain$results$ROC,
	'\n')


Predd <- predict(gbmTrain, newdata=testing, type="prob")
colAUC(Predd, testing, plot=T)

## COMPUTE A CONFUSION MATRIX
confusionMatrix(data=gbmTest, testing$species)




#### Redundancy analysis (RDA) (C)
set.seed(34)
#training
rdaTrain <- train(species~.,
			data=training,
			method="rda",
			tuneGrid=data.frame(.lambda=1, .gamma=0),
			trControl=trainControl(method="repeatedcv",	number=10, repeats=5,
					classProbs=T,
					summaryFunction=twoClassSummary		## ROC curve
					),
			preProc=c("center","scale"),
			metric="ROC"
			)
rdaTrain

plot(rdaTrain)

## RMSE + R2
rdaTest <- predict(rdaTrain, newdata=testing)
cat('\nComputeRMSE','\nRMSE=',
	sqrt((sum((testing$species-rdaTest)^2))/nrow(testing)),
	'\n',
	'\nCompute Rsquared', '\nR2=',
	rdaTrain$results$Rsquared,
	'\n')

#### For binary classification (only)
## AUC
rdaTest <- predict(rdaTrain, newdata=testing)
cat('\nComputeROC........','ROC=',
	rdaTrain$results$ROC,
	'\n')


Predd <- predict(rdaTrain, newdata=testing, type="prob")
colAUC(Predd, testing, plot=T)

## COMPUTE A CONFUSION MATRIX
confusionMatrix(data=rdaTest, testing$species)





#### Linear regression (START)
#### Partial least squares (D)
set.seed(34)
#training
plsTrain <- train(species~.,
			data=training,
			method="pls",
			trControl=trainControl(method="repeatedcv",	number=10, repeats=5
					#classProbs=T,
					#summaryFunction=twoClassSummary		## ROC curve
					),
			preProc=c("center","scale")
			#metric="ROC"
			)
plsTrain

plot(plsTrain)

## RMSE + R2
plsTest <- predict(plsTrain, newdata=testing)
cat('\nComputeRMSE','\nRMSE=',
	sqrt((sum((testing$species-plsTest)^2))/nrow(testing)),
	'\n',
	'\nCompute Rsquared', '\nR2=',
	plsTrain$results$Rsquared,
	'\n')

#### For binary classification (only)
Predd <- predict(plsTrain, newdata=testing, type="prob")
colAUC(Predd, testing, plot=T)

## COMPUTE A CONFUSION MATRIX
confusionMatrix(data=plsTest, testing$species)



#### Multinomial logistic regression (MLR) (C)
set.seed(34)
#training
multiTrain <- train(species~.,
			data=training,
			method="multinom",
			trControl=trainControl(method="repeatedcv",	number=10, repeats=5
					#classProbs=T,
					#summaryFunction=twoClassSummary		## ROC curve
					),
			preProc=c("center","scale")
			#metric="ROC"
			)
multiTrain

plot(multiTrain)
## RMSE + R2
multiTest <- predict(multiTrain, newdata=testing)
cat('\nComputeRMSE','\nRMSE=',
	sqrt((sum((testing$species-multiTest)^2))/nrow(testing)),
	'\n',
	'\nCompute Rsquared', '\nR2=',
	multiTrain$results$Rsquared,
	'\n')

#### For binary classification (only)
Predd <- predict(multiTrain, newdata=testing, type="prob")
colAUC(Predd, testing, plot=T)

## COMPUTE A CONFUSION MATRIX
confusionMatrix(data=multiTest, testing$species)





#### GLM (D)
set.seed(34)
#training
glmTrain <- train(species~.,
			data=training,
			method="glm",
			trControl=trainControl(method="repeatedcv",	number=10, repeats=5
					#classProbs=T,
					#summaryFunction=twoClassSummary		## ROC curve
					),
			preProc=c("center","scale")
			#metric="ROC"
			)
glmTrain

plot(glmTrain)

## RMSE + R2
glmTest <- predict(glmTrain, newdata=testing)
cat('\nComputeRMSE','\nRMSE=',
	sqrt((sum((testing$species-glmTest)^2))/nrow(testing)),
	'\n',
	'\nCompute Rsquared', '\nR2=',
	glmTrain$results$Rsquared,
	'\n')

#### For binary classification (only)
Predd <- predict(glmTrain, newdata=testing, type="prob")
colAUC(Predd, testing, plot=T)

## COMPUTE A CONFUSION MATRIX
confusionMatrix(data=glmTest, testing$species)




#### GLM boostrapping (D)
set.seed(21)
#training
gbTrain <- train(species~.,
			data=training,
			method="glmboost",
			tuneGrid=data.frame(.mstop=50, .prune='0'),
			trControl=trainControl(method="repeatedcv",	number=10, repeats=5
					#classProbs=T,
					#summaryFunction=twoClassSummary		## ROC curve
					),
			preProc=c("center","scale")
			#metric="ROC"
			)
gbTrain

plot(gbTrain)

## RMSE + R2
gbTest <- predict(gbTrain, newdata=testing)
cat('\nComputeRMSE','\nRMSE=',
	sqrt((sum((testing$species-gbTest)^2))/nrow(testing)),
	'\n',
	'\nCompute Rsquared', '\nR2=',
	gbTrain$results$Rsquared,
	'\n')

#### For binary classification (only)
Predd <- predict(gbTrain, newdata=testing, type="prob")
colAUC(Predd, testing, plot=T)

## COMPUTE A CONFUSION MATRIX
confusionMatrix(data=gbTest, testing$species)




#### Bayesian GLM (D)
set.seed(21)
#training
bglmTrain <- train(species~.,
			data=training,
			method="bayesglm",
			trControl=trainControl(method="repeatedcv",	number=10, repeats=5
					#classProbs=T,
					#summaryFunction=twoClassSummary		## ROC curve
					),
			preProc=c("center","scale")
			#metric="ROC"
			)
bglmTrain

plot(bglmTrain)

## RMSE + R2
bglmTest <- predict(bglmTrain, newdata=testing)
cat('\nComputeRMSE','\nRMSE=',
	sqrt((sum((testing$species-bglmTest)^2))/nrow(testing)),
	'\n',
	'\nCompute Rsquared', '\nR2=',
	bglmTrain$results$Rsquared,
	'\n')

#### For binary classification (only)
Predd <- predict(bglmTrain, newdata=testing, type="prob")
colAUC(Predd, testing, plot=T)

## COMPUTE A CONFUSION MATRIX (very important cf. 04-7 slide 32 of ch4 statlearning)
confusionMatrix(data=bglmTest, testing$species)



## PREDICTION OF NEW SAMPLES (BATCH COMPARISON)
models <- list(svm=svmTrain, nnet=nnetTrain, RF=rfTrain, gbm=gbmTrain)
testPredd <- predict(models, newdata=testing[,-dim(testing)[2]])
lapply(testPredd, function(x) x[100:110])
testing[100:110,]

## HOW DO THESE MODELS COMPARE TO THEIR RESAMPLING RESULTS
comp.mods <- resamples(list(rda=rdaTrain, svm=svmTrain,nnet=nnetTrain))
summary(comp.mods)
# plot confusion rates
x11(); xyplot(comp.mods, what="BlandAltman")
## assess the difference in the AUC
diffs <- diff(comp.mods)
summary(diffs)



setwd("C:/Dropbox/Workshop2013/Work/R/ANN/")
lsos(pat="Train.*")
save(list=ls(pattern="Train.*"),file="caret.Rdata")


=========================================================================
#	SECTION --								#
=========================================================================
## ENSEMBLE LEARNING

library(foreach)
library(doSNOW)
getDoParWorkers()
getDoParName()
getDoParVersion()

cl <- makeCluster(4)
registerDoSNOW(cl)

ctl=expand.grid(.C=2.2, .sigma=0.27)	## svmRadial
ctl=expand.grid(.size=1, .decay=5e-6)	## nnet
ctl=expand.grid(.mtry=(17))	## Random forest
ctl=expand.grid(.mstop=50, .prune='0')	## GLM boostrap
ctl=expand.grid(.interaction.depth=c(1,5,9), .n.trees=(1:15)*100,.shrinkage=0.1)

## BEGIN RUN
set.seed(21)
system.time(
svmbag <- bagging(training,testing,m=10,ite=5000,meth="nnet",gridZ=ctl)
)
## END RUN
sqrt((sum((testing$species - svmbag)^2))/nrow(testing))		## compute RMSE

stopCluster(cl)		## close cluster only after finishing w all models

## COMPUTE RMSE
ensemble.pred <- (gbbag+rfbag)/2
ensemble.pred <- (gbbag*2+rfbag)/3
ensemble.pred <- (gbbag+rfbag*2)/3
sqrt((sum((testing$species - ensemble.pred)^2))/nrow(testing))


## plot RMSE for the different learners
## copy/paste content of *.xls file
## source file can be found here : C:\Users\Sleiman Bassim\Downloads\Learners
x <- read.table("clipboard", sep="\t", header=T);x
library(lattice)

## plot RMSE versus different learners
xyplot(SumRMSE ~ Learners | weight + Bagging, data=x, type="o", col="black",pch=20, cex=1)

## plot the duration of the run with different clusters (4 and 7 of sNOW)
xyplot(Time ~ Iterations | Clusters + Learner, data=x, type="o", col="blue",pch=20, cex=1)

=========================================================================
#	SECTION --								#
=========================================================================
## using standalone packages

require(randomForest)
rff <- randomForest(species~., data=training, ntree=500)
rd.pred <- predict(rff, newdata=testing)

bin <- class.ind(training[,dim(training)[2]])
irisANN <- nnet(training[,-dim(training)[2]], bin, size=1, softmax=T)
nnetTest <- predict(irisANN, newdata=testing)

=========================================================================
#								SECTION --								#
=========================================================================
## COMPARE LEARNERS (WITHIN EACH MODEL GROUP)
set.seed(21)
ctrl <- trainControl(method="repeatedcv", number=10, repeats=5)

## Try classifiers
svmRad <- train(species~., data=NN, tuneLength=10, method="svmRadial", trControl=ctrl, preProc=c("center","scale"))
svmPoly <- train(species~., data=NN, tuneLength=10, method="svmPoly", trControl=ctrl, preProc=c("center","scale"))
svmLin <- train(species~., data=NN, tuneLength=10, method="svmLinear", trControl=ctrl, preProc=c("center","scale"))
svmWei <- train(species~., data=NN, tuneLength=10, method="svmRadialWeights", trControl=ctrl, preProc=c("center","scale"))

## Compare
compare.svm <- resamples(list(Radial=svmRad, Poly=svmPoly, Linear=svmLin))
summary(compare.svm)
#bwplot(compare.svm metric="Accuracy")
#densityplot(compare.svm, metric="Accuracy")
xyplot(compare.svm, what="BlandAltman")

=========================================================================
#	SECTION --								#
=========================================================================

## adaboost package for boosting
require(ada)
data(Smarket, package="ISLR")
x = Smarket
n = dim(x)[1]
Tr.index <- sample(1:n, floor(0.6*n), FALSE)
training <- x[Tr.index, ]
testing <- x[-Tr.index, ]

gdis <- ada(Direction ~., data=training, iter=20, nu=1, type="discrete")
gdis
gdis <- addtest(gdis, testing[,-dim(testing)[2]], testing[,dim(testing)[2]])
plot(gdis, T, T)
varplot(gdis)
pairs(gdis, training[,-dim(training)[2]],maxvar=2)
gTest <- predict(gdis, newdata=testing[,-dim(testing)[2]])


=========================================================================
#	SECTION --								#
=========================================================================
## Removing features w more than 0.7 correlation

library(corrplot)
library(caret)
data(Smarket, package="ISLR")
x=data(iris)
x=Smarket
data.scale <- scale(x[,-dim(x)[2]], center=TRUE, scale=T)
head(data.scale)
data.cor <- cor(data.scale)	## compute the correlation matrix
head(data.cor)
corrplot(data.cor, order="hclust")
high.cor <- findCorrelation(data.cor, 0.70)	## find highly correlated features
data.scale2 <- data.scale[,-high.cor]		## remove highly correlated features
data.cor <- cor(data.scale2)
corrplot(data.cor, order="hclust")







====================
    stat LEARNING
====================

## source Video 0.4.1 Classification in R ch4 from StatLearning
data(Smarket, package="ISLR") ## binary categorical
pairs(Smarket, col=Smarket$Direction)
cor(Smarket[,-9])
glm.fit <- glm(Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume, data = Smarket, family = binomial)
summary(glm.fit)
glm.probs <- predict(glm.fit, type="response")
glm.pred <- ifelse(glm.probs>0.5, "Up","Down")
attach(Smarket)
table(glm.pred, Direction)
mean(glm.pred==Direction)
par(mfrow = c(2,2))
plot(glm.fit)
## prediction for a binomial classification (logistic regression)

train <- Year<2005
str(train)
glm.fit <- glm(Direction~Lag1+Lag2, data = Smarket, family = binomial, subset=train)
glm.probs <- predict(glm.fit, newdata=Smarket[!train,], type = "response")
glm.pred <- ifelse(glm.probs>0.5, "Up","Down")
table(glm.pred, Direction[!train])
mean(glm.pred==Direction[!train])
## make a training and testing set

require(MASS)
lda.fit <- lda(Direction~Lag1+Lag2, data=Smarket, subset=Year<2005)
lda.fit
plot(lda.fit)
Smarket.2005 <- subset(Smarket, Year==2005)
lda.pred <- predict(lda.fit, Smarket.2005)
head(data.frame(lda.pred))
table(lda.pred$class, Smarket.2005$Direction)
mean(lda.pred$class==Smarket.2005$Direction)
## Linear Discriminat analysis LDA

require(class)
Xlag <- cbind(Lag1, Lag2)
train <- Year<2005
knn.pred <- knn(Xlag[train,], Xlag[!train,], Direction[train], k=3)
table(knn.pred, Direction[!train])
mean(knn.pred==Direction[!train])
detach(Smarket)
## K-nearest neighbor

library(ISLR)
set.seed(123)
train <- sample(1:392,196)
lm.fit <- lm(mpg~horsepower, data = Auto, subset = train)
attach(Auto)
mean((mpg-predict(lm.fit, Auto))[-train]^2)
lm.fit2 <- lm(mpg~poly(horsepower,2), data = Auto, subset=train)
mean((mpg-predict(lm.fit2, Auto))[-train]^2)
lm.fit3 <- lm(mpg~poly(horsepower,3), data = Auto, subset=train)
mean((mpg-predict(lm.fit3, Auto))[-train]^2)
##the validation set approach

require(boot)
glm.fit <- glm(mpg~horsepower, data=Auto)	## i used glm to do a linear regression & is used w/ cv.glm
cv.err <- cv.glm(Auto, glm.fit)
cv.err$delta	## error
# example 1 of LOOCV

cv.err <- rep(1:5)
for(i in 1:5){
glm.fit <- glm(mpg~poly(horsepower, i), data = Auto)
cv.err[i] <- cv.glm(Auto, glm.fit)$delta[1]
}
cv.err
## Leave one out cross validation (LOOCV)

cv.err <- rep(1:10)
for(i in 1:10){
glm.fit <- glm(mpg~poly(horsepower, i), data = Auto)
cv.err[i] <- cv.glm(Auto, glm.fit, K=10)$delta[1]
}
cv.err
## K-fold cross validation

require(boot)
boot.fn <- function(data, index){
return(coef(lm(mpg~horsepower, data = Auto, subset=index)))}
boot.fn(Auto, 1:392)
boot.fn(Auto, sample(392,392, replace = T))
fit.boot <- boot(Auto, boot.fn, 1000)
plot(fit.boot)
## Bootstrap (linear regression)

boot.fn <- function(data, index)
    return(coef(lm(mpg~horsepower+I(horsepower^2), data = Auto, subset=index)))
set.seed(1)
boot(Auto, boot.fn, 1000)
summary(lm(mpg~horsepower+I(horsepower^2),data = Auto))$coefficients
## Bootstrap (polynomial)

require(ISLR)
require(boot)
cv.err <- rep(1:10)
degree <- 1:10
for(i in degree){
glm.fit <- glm(Volume~poly(Lag1, i), data = Smarket)
cv.err[i] <- cv.glm(Smarket, glm.fit, K=10)$delta}
cv.err
plot(degree,cv.err, type="b", ylab="Error Rate", xlab="Polynomial degree",col="chocolate")


rm(list=ls())
setwd("c:/Users/Sleiman Bassim/Downloads/")
load("5.R.Rdata", .GlobalEnv)
lsos(pat="")
summary(Xy)
attach(Xy)
summary(lm(y~X1+X2))
matplot(Xy, type="l")
boot.fn <- function(data, index) return(coef(lm(y~X1+X2,data=data, subset=index)))
boot.fn(Xy, 1:dim(Xy)[1])
require(boot)
boot(Xy, boot.fn, 100)
detach(Xy)
## Bootstrap standard

set.seed(123)
boot.fn <- function(data) return(coef(lm(y~X1+X2,data=data)))
tsboot(Xy, boot.fn, 1000, sim = "fixed", l=100)
## bootstrap blocks


require(boot)
set.seed(123)
x1 <- rnorm(100)
var(x1);mean(x1)
x3 <- rnorm(100)
x2 <- rnorm(100)
y <- x1-2*x2^2+x3^3+rnorm(100)
par(mfrow = c(2,2))
plot(y~x1+x2+x3)
norm.data <- as.data.frame(cbind(x1,x2,x3,y))
summary(norm.data)
fit.glm <- glm(y~., data=norm.data)
cv.err <- cv.glm(norm.data, fit.glm, K=5)$delta[1]
cv.err			## return the standard erro estimates of the linear model accuracy
boot.fn <- function(data, index)return(coef(glm(y~., data=data, subset = index)))
boot.fn(norm.data, 1:100)	## return the coeficients of the linear model
boot(norm.data, boot.fn, 1000)	## bootstrap; return std errors only
summary(glm(y~., data = norm.data))	## return coeficents and standard errors
cv.err <- 1:10
index <- 1:10
for(i in index){
glm.fit <- glm(y~poly(x1+x2+x3, i),data = norm.data)
cv.err[i] <- cv.glm(norm.data, glm.fit, K=10)$delta[1]
}
plot(index, cv.err, xlab = "Polynomial degree", type = "b")
## Simulated data set for cross validation and bootstrap

require(ISLR)
names(Hitters)
dim(Hitters)
with(Hitters,sum(is.na(Salary)))
Hitters <- na.omit(Hitters)
dim(Hitters)
## remove NAs
require(leaps)
regfit.full <- regsubsets(Salary~., Hitters, nvmax = 19)
## best subsets on all 19 variables
reg.summary <- summary(regfit.full)
names(reg.summary)
plot(reg.summary$bic, type="b")
points(which.min(reg.summary$bic), min(reg.summary$bic),col="red",pch=20)
?plot.regsubsets
plot(regfit.full, scale="bic")
coef(regfit.full, 6)	## look at the coefficient estimates of the six-variable selected model
## best subset selection

regfit.fwd <- regsubsets(Salary~., Hitters, nvmax = 19, method = "forward")
reg.summary <- summary(regfit.fwd)
names(reg.summary)
plot(reg.summary$bic, type="b")
points(which.min(reg.summary$bic), min(reg.summary$bic),col="red",pch=20)
plot(regfit.fwd, scale="bic")
## Forward stepwise selection


regfit.bwd <- regsubsets(Salary~., Hitters, nvmax = 19, method = "backward")
reg.summary <- summary(regfit.bwd)
names(reg.summary)
plot(reg.summary$bic, type="b")
points(which.min(reg.summary$bic), min(reg.summary$bic),col="red",pch=20)
plot(regfit.fwd, scale="bic")
## Backward stepwise selection
coef(regfit.fwd, 7)
coef(regfit.full, 7)
coef(regfit.bwd, 7)
## the 7th model is different for all 3 subset selection methods, even if the 6ths are identical

set.seed(1)
x <- dim(Hitters)[1]
train <- sample(x, 2*x/3,rep=T)
regfit.best <- regsubsets(Salary~., Hitters[train,], nvmax = 19)
test.mat <- model.matrix(Salary~.,data = Hitters[!train,])
## cross validation
require(glmnet)
require(ISLR)
sum(is.na(Hitters$Salary))
Hitters <- na.omit(Hitters)
x <- model.matrix(Salary~., Hitters)[,-1]
y <- Hitters$Salary
## Prepare data for glmnet()
grid <- 10^seq(10,-2,length=100)
ridge.mod <- glmnet(x, y, alpha=0, lambda=grid, standardize = T)	## alpha=1 lasso =0 ridge regression
dim(coef(ridge.mod))
## shrinkage of 20 predictors over 100 lambda
names(ridge.mod)
ridge.mod$lambda[50]
coef(ridge.mod)[,50]
sqrt(sum(coef(ridge.mod)[-1,50]^2))
sqrt(sum(coef(ridge.mod)[-1,60]^2))
## compute l2 norm of the coefficients
predict(ridge.mod, s=50, type = "coefficients")[1:20,]
## predict new coefficients for lambda=50
set.seed(2)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)
y.test <- y[test]
## prepare datasets to compute test MSE of the ridge regression
ridge.mod <- glmnet(x[train,], y[train], alpha = 0, lambda=grid, thresh = 1e-12)
ridge.pred <- predict(ridge.mod, s=4, newx=x[test,])
mean((ridge.pred-y.test)^2)
ridge.pred <- predict(ridge.mod, s=0, newx=x[test,], exact=T)	## to predict with a lambda of 0
mean((ridge.pred-y.test)^2)
## test different lambdas
lm(y~x, subset = train)
predict(ridge.mod, s=0, exact=T, type="coefficients")[1:20,]
## lm and lambda=0 are the same, but lm() provides more useful outputs
set.seed(3)
cv.out <- cv.glmnet(x[train,], y[train],alpha=0, nfolds = 10)
plot(cv.out)
## cross validation on ridge regression
best.lambda <- cv.out$lambda.min
best.lambda
ridge.pred <- predict(ridge.mod, s=best.lambda, newx=x[test,])
mean((ridge.pred-y.test)^2)
## choose the best lambda by cross validation
out <- glmnet(x, y, alpha = 0)
predict(out, type="coefficients", s=best.lambda)[1:20,]
## refit the model using all the dataset
## Ridge regression (shrinkage)

lasso.mod <- glmnet(x[train,], y[train], alpha = 1, lambda=grid)
plot(lasso.mod)
set.seed(3)
cv.out <- cv.glmnet(x[train,], y[train],alpha=1)
plot(cv.out)
bestlambda <- cv.out$lambda.min
lasso.pred <- predict(lasso.mod, s=bestlambda, newx=x[test,])
mean((lasso.pred-y[test])^2)
out <- glmnet(x, y, alpha=1)
lasso.coef <- predict(out, s=bestlambda, type="coefficients")[1:20,]
lasso.coef
lasso.coef[lasso.coef!=0]
## lasso (shrinkage)

x <- rnorm(100)
y <- 2+x^2
fit <- lm(y~x)
plot(y~x, col="chocolate",pch=20)
abline(fit,col="forestgreen",lwd=1)
pred <- predict(fit)
segments(x,y,x,pred, lty="dashed")
points(x,pred,pch=4,cex=.5)
## simulated linear model data


require(pls)
set.seed(4)
x <- model.matrix(Salary~., Hitters)[,-1]
y <- Hitters$Salary
pcr.fit <- pcr(y~x, scale=T, validation="CV")
summary(pcr.fit)
## calculating the RMSE for each component
350.9^2
## calculating the MSE
validationplot(pcr.fit, val.type = "MSEP")
## ploting the MSE for each component
train <- sample(1:nrow(Hitters), nrow(Hitters)/2)
test <- (-train)
pcr.fit <- pcr(y~x, subset=train, scale=T, validation="CV")
validationplot(pcr.fit, val.type = "MSEP")
## MSE for the training data
pcr.pred <- predict(pcr.fit, x[test,], ncomp=5)	## 5 is the lowest component MSE
mean((pcr.pred-y[test])^2)
## MSE for the test set
pcr.fit <- pcr(y~x, scale=TRUE, ncomp=5)
summary(pcr.fit)
## fit the moedl on all the data set
## Principal component regression


require(pls)
set.seed(5)
x <- model.matrix(Salary~., Hitters)[,-1]
y <- Hitters$Salary
train <- sample(1:nrow(Hitters), nrow(Hitters)/2)
test <- (-train)
pls.fit <- plsr(y~x, subset=train, scale=T, validation="CV")
summary(pls.fit)
validationplot(pls.fit, val.type = "MSEP")
pls.pred <- predict(pls.fit, x[test,], ncomp=2)
mean((pls.pred-y[test])^2)
pls.fit <- plsr(y~x, ncomp=2)
summary(pls.fit)
## Partial least squares









