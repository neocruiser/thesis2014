
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
xyplot(Profile$results$ROC ~ Profile$results$Variables, type=c("p","g"),auto.key=T)
Profile$fit
head(Profile$resample)
plot(Profile)


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
x = Smarket				## original dataset
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

## COMPUTE A CONFUSION MATRIX
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


















