#########################################################################
#								SECTION -1-								#
#########################################################################
## source: http://tinyurl.com/owvbekw
## CLASSIFICATION

require(nnet)
setwd("C:\\Dropbox\\Workshop2013\\Work\\R\\ANN\\")
seeds <- read.table('seeds.txt', sep="\t")

seedstrain<- sample(1:210,147)
seedstest <- setdiff(1:210,seedstrain)
ideal <- class.ind(seeds$V8)
seedsANN = nnet(seeds[seedstrain,-8], ideal[seedstrain,], size=10, softmax=TRUE)
predict(seedsANN, seeds[seedstrain,-8], type="class")
table(predict(seedsANN, seeds[seedstest,-8], type="class"), seeds[seedstest,]$V8)

#########################################################################
#								SECTION -2-								#
#########################################################################
## BASIC CLASSIFICATION IN R: NEURAL NETWORKS AND SUPPORT VECTOR MACHINES
## source : http://tinyurl.com/k3gbsuo

require(kernlab)
data(package="kernlab")

## LOAD IRIS DATA
## ORIGINAL SETTINGS SIMILAR TO THE SOURCE
irisdata <- read.csv("http://www.heatonresearch.com/dload/data/iris.csv",head=TRUE,sep=",")
summary(irisdata)
## TRAINING AND VALIDATION
## create vectors
irisTrain <- sample(1:150,100)
irisVal <- setdiff(1:150, irisTrain)
## USING SUPPORT VECTOR MACHINE (SVM) for training
## create a radial basis function(RBF)
rbf <- rbfdot(sigma=0.1)
## train the SVM
irisSVM <- ksvm(species ~ ., data=irisdata[irisTrain,], type="C-bsvc", kernel=rbf, C=10, prob.model=T)
fitted(irisSVM)		## the fitted values
## validate with a prediction, while removing the 5th col and predicting the species
predict(irisSVM, irisdata[irisVal,-5], type="probabilities")

## FROM THE MANUAL OF KERNLAB
require(kernlab)
irisdata <- read.csv("http://www.heatonresearch.com/dload/data/iris.csv",head=TRUE,sep=",")
irisSample <- sample(1:dim(irisdata)[1],5)
irisTrain <- irisdata[-irisSample,]
irisTest <- irisdata[irisSample,]
gene <- ksvm(species ~., data=irisTrain, kernel="rbfdot", kpar="automatic", C=60, cross=3, prob.model=T)
gene
predict(gene,irisTest)
predict(gene, irisTest, type="probabilities")
## TOY CLASSIFICATION WITH BINARY CLASSES
x <- rbind(matrix(rnorm(120), ,2), matrix(rnorm(120, mean=3),,2))
y<- matrix(c(rep(1,60), rep(-1,60)))
svp <- ksvm(x,y,type="C-svc")
plot(svp, data=x)
plot(ran)
plot(ran, cex=c[ranked[,3]]/40)

## RANKING
require(kernlab)
data(spirals)
head(spirals);tail(spirals)
ran <- spirals[rowSums(abs(spirals) < 0.55) == 2,]
ranked <- ranking(ran, 54, kernel="rbfdot", kpar=list(sigma=100), edgegraph=T)
ranked[54,2] <- max(ranked[-54,2]) ## remplacer un NA
c=1:86
op <- par(mfrow=c(1,2), pty="s")
plot(ran)
plot(ran,cex=c[ranked[,3]]/40)

## USING AN ANN
## SIMILAR TO THE SOURCE
require(nnet)
irisdata <- read.csv("http://www.heatonresearch.com/dload/data/iris.csv",head=TRUE,sep=",")
summary(irisdata)
## TRAINING AND VALIDATION
## create vectors
irisTrain <- sample(1:150,100)
irisVal <- setdiff(1:150, irisTrain)
## NORMALIZE ONE OF N (0 TO 1)
ideal <- class.ind(irisdata$species)
## TRAIN A NN
irisANN <- nnet(irisdata[irisTrain,-5], ideal[irisTrain,], size=10, softmax=T)
## VALIDATE
predict(irisANN, irisdata[irisVal,-5],type="class")

#########################################################################
#								SECTION -3-								#
#########################################################################
 ## source: http://tinyurl.com/n57pn7a
require(neuralnet)
## binarize the categorical output
nnet_iristrain <- irisdata
nnet_iristrain <- cbind(nnet_iristrain, irisdata$species == 'Iris-setosa', irisdata$species == 'Iris-versicolor', irisdata$species == 'Iris-virginica')
names(nnet_iristrain) <- c(colnames(nnet_iristrain[1:5]), "setosa","versicolor","virginica")
nnet_iristrain
nn <- neuralnet(setosa+versicolor+virginica ~ sepal_l+sepal_w+petal_l+petal_w, data=nnet_iristrain, hidden=c(3))
plot(nn)

mypredict <- compute(nn, iristest[-5])$net.result


#########################################################################
#								SECTION -4-								#
#########################################################################
## K CROSS VALIDATION or adjR2 or BIC model selection (before training)
## http://tinyurl.com/pyye5x6

irisdata <- read.csv("http://www.heatonresearch.com/dload/data/iris.csv",head=TRUE,sep=",")
require(kernlab)
require(leaps)
set.seed(28)

## FIRST APPROACH
## MODEL SELECTION ACCORDING TO ADJR2 OR BIC
model <- regsubsets(species~., data=irisdata)
sum.model <- summary(model)
sum.model
names(sum.model)
## PLOT THE MAXIMUM ADJUSTE R2
## According to adjR2
plot(sum.model$adjr2, xlab="Number of variables", ylab="adj-R2", type='l')
x <- which.max(sum.model$adjr2)
points(x, sum.model$adjr2[x], col="red", cex=2, pch=20)
## Residual SUM OF SQUARES
plot(sum.model$rss, xlab="Number of variables", ylab="BIC", type='l')
## PLOT USING THE REGSUBSETS FUNCTION
names(sum.model)
plot(model, scale="r2")
plot(model, scale="adjr2")
plot(model, scale="Cp")
plot(model, scale="bic") ## according to BIC
## COEFFICIENTS
coef(model, 4) ## best model '4'

## SECOND APPROACH
## MODEL SELECTION ACCORDING VALIDATION AND CROSSVALIDATION
irisdata01 <- data.frame(irisdata[,-5],species=c(rep(0,50), rep(1,50), rep(2,50)))
set.seed(134)
train <- sample(c(TRUE, FALSE), nrow(irisdata01), rep=T)
test <- (!train)
# fit the model using the training set
model.best <- regsubsets(species~., data=irisdata01[train,]) 	# nvmax= ?
sum.model.best <- summary(model.best)
sum.model.best
## Predict on coefficients of the test set
## FUNCTION
predict.regsubsets <- function(object, newdata, id, ...){
	form <- as.formula(object$call[[2]])
	mat <- model.matrix(form, newdata)
	coefi <- coef(object, id=id)
	xvars <- names(coefi)
	mat[,xvars] %*%coefi
	}

## CROSSVALIDATION FOR PREDICTING THE NUMBER OF VARIABLES ARE BEST FOR GENERALIZATION
k=10
set.seed(123)
folds=sample(1:k, nrow(irisdata01), replace=T)
## create a matrix for CV
cv.errors <- matrix(NA, k, 4, dimnames=list(NULL, paste(1:4)))
for(j in 1:k){
	best.fit <- regsubsets(species~., data=irisdata01[folds!=j,])
	for(i in 1:4){
		pred <- predict(best.fit, irisdata01[folds==j,], id=i)
	cv.errors[j,i] = mean((irisdata01$species[folds==j]-pred)^2)
	}
}

mean.cv.errors <- apply(cv.errors, 2, mean)
mean.cv.errors
which.min(mean.cv.errors)
plot(mean.cv.errors, type="b")
## extract the coefficients on all data set
reg.best <- regsubsets(species~.,data=irisdata01)
coef(reg.best, 3)
summary(reg.best)


###### USING SUPPORT VECTOR MACHINE (SVM) FOR TRAINING THEN PREDICTION
## create vectors
require(kernlab)
train <- sample(1:dim(Smarket)[1], 1000)
test <- setdiff(1:dim(Smarket)[1], train)
## create a radial basis function(RBF)
rbf <- rbfdot(sigma=0.1)
## train the SVM
irisSVM <- ksvm(Year ~ ., data=Smarket[train,-2], type="C-bsvc", kernel="rbfdot", C=10, prob.model=T)
fitted(irisSVM)		## the fitted values
## validate with a prediction on the test data (qualification)
predict(irisSVM, Smarket[test,-c(1,2)], type="probabilities")

##### USING AN ANN FOR TRAINING AND PREDICTION
require(nnet)
train <- sample(1:dim(Smarket)[1], 1000);length(train)
test <- setdiff(1:dim(Smarket)[1], train);length(test)
## NORMALIZE ONE OF N (0 TO 1)
ideal <- class.ind(Smarket$Year);dim(ideal)
AA <- data.frame(Smarket[,-1], ideal)
## TRAIN A NN
irisANN <- nnet(Year~., data=Smarket[train,], size=10, softmax=F, decay=0, maxit=2000, trace=T, linout=T)
## VALIDATE
predict(irisANN, Smarket[test,-1],type="class")


#########################################################################
#								SECTION -2-								#
#########################################################################
## MODEL SELECTION (before training) Feature selection
## http://tinyurl.com/pyye5x6

require(leaps)


## FIRST APPROACH
##### MODEL SELECTION ACCORDING TO ADJR2 OR BIC
model <- regsubsets(Direction~., data=Smarket, really.big=F, nvmax=9)
sum.model <- summary(model)
sum.model
names(sum.model)
## PLOT THE MAXIMUM ADJUSTE R2
## According to adjR2
plot(sum.model$adjr2, xlab="Number of variables", ylab="adj-R2", type='l')
x <- which.max(sum.model$adjr2)
points(x, sum.model$adjr2[x], col="red", cex=2, pch=20)
## Residual SUM OF SQUARES
plot(sum.model$rss, xlab="Number of variables", ylab="BIC", type='l')
## PLOT USING THE REGSUBSETS FUNCTION
plot(model, scale="adjr2")
plot(model, scale="bic") ## according to BIC

## SECOND APPROACH
###### CROSSVALIDATION FOR PREDICTING THE NUMBER OF VARIABLES ARE BEST FOR GENERALIZATION
set.seed(134)
train <- sample(c(TRUE, FALSE), nrow(Smarket), rep=T)
test <- (!train)
# fit the model using the training set
model.best <- regsubsets(Direction~., data=Smarket[train,], nvmax=9) 	# nvmax= ?
sum.model.best <- summary(model.best)
sum.model.best
## Predict on coefficients of the test set
## FUNCTION
predict.regsubsets <- function(object, newdata, id){
	form <- as.formula(object$call[[2]])
	mat <- model.matrix(form, newdata)
	coefi <- coef(object, id=id)
	xvars <- names(coefi)
	mat[,xvars] %*%coefi
	}

## MODEL SELECTION ACCORDING VALIDATION AND CROSSVALIDATION
k=10
set.seed(123)
folds=sample(1:k, nrow(Smarket), replace=T)
## create a matrix for CV
cv.errors <- matrix(NA, k, 8, dimnames=list(NULL, paste(1:8)))
for(j in 1:k){
	best.fit <- regsubsets(Direction~., data=Smarket[folds!=j,], nvmax=8) # train # change ...~
	for(i in 1:8){
		pred <- predict(best.fit, Smarket[folds==j,], id=i) # test
	cv.errors[j,i] = mean((Smarket$Direction[folds==j]-pred)^2)	## change $...
	}
}

mean.cv.errors <- apply(cv.errors, 2, mean)
mean.cv.errors
which.min(mean.cv.errors)
plot(mean.cv.errors, type="b")
## extract the coefficients using all the data set
reg.best <- regsubsets(Direction~.,data=Smarket)
coef(reg.best, 7)
summary(reg.best)




#########################################################################
#								SECTION --								#
#########################################################################
## ENSEMBLE LEARNING

library(foreach)
library(caret)
library(doSNOW)

## Bootstrap aggregation (Bagging) + 1 classifier = RMSE
bagging <- function(training, testing, m=10, ite=2, meth="svmRadial", tune=10){

	Predd <- foreach(i=1:ite,.combine=cbind,.export=c("train","trainControl")) %dopar% {
	
		bagging.index <- sample(1:dim(training)[1], size=floor((nrow(training)/m)))	## vector list
					
		Train.me <- train(species~., 
				data=training[bagging.index,],
				method=meth,
				trControl=trainControl(method="repeatedcv", number=10, repeats=3
					#classProbs=T,
					#summaryFunction=twoClassSummary
					),
				tuneLength=tune,
				preProc=c("center","scale")
				#metric="ROC"
				)
			
	predict(Train.me, newdata=testing)
}
	rowMeans(Predd)
}

data(Smarket)
x = Smarket				## original dataset
NN <- class.ind(x[,dim(x)[2]])
NN <- data.frame(x[,-dim(x)[2]], species=NN[,1])		## all dataset binary
y = NN[,dim(NN)[2]]
## CREATE CROSS VALIDATION K FOLDS FOR TRAINING
trainList <- createDataPartition(y, p = 0.8, list = FALSE)
training <- NN[trainList,]
testing <- NN[-trainList,]

## Begin run
cl <- makeCluster(2)
registerDoSNOW(cl)
svmPoly <- bagging(training,testing,m=10,ite=50,meth="svmPoly",tune=10)
stopCluster(cl)
## End run
ensemble.pred <- (svmPred+nnetPred)/2
sqrt((sum((testing$species - ensemble.pred)^2))/nrow(testing))



## Slow run
cl <- makeCluster(2)
registerDoSNOW(cl)
set.seed(6) ; ite <- 50 ; m <- 10
Predd <- foreach(i=1:ite, .combine=cbind, .export=c("train","trainControl")) %dopar% {
	
	bagging.index <- sample(1:dim(training)[1], size=floor((nrow(training)/m)))
	#bagging.pos <- 1:nrow(training) %in% bagging.index
	svmTrain <- train(species~., 
			data=training[bagging.index,],
			method="svmRadial",
			trControl=trainControl(method="repeatedcv", number=10, repeats=3
					#classProbs=T,
					#summaryFunction=twoClassSummary
					),
			preProc=c("center","scale"),
			tuneLength=10
			#metric="ROC"
			)
	predict(svmTrain, newdata=testing)
}
stopCluster(cl)
svmPred <- rowMeans(Predd)
sqrt((sum((testing$species - svmPred)^2))/nrow(testing))
















