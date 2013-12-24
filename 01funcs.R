## Functions for R
## udpated on Friday, 12/2013
## load with C-c C-l in emacs before sessions


## EBDBN
## FUNCTION THAT TAKE THE ZSCORES $Z AND CREATE A BINARY MATRIX WITH 1 AND -1 FOR UP AND DOWN REGULATED DET FOR 0.99 THRESHOLD
net2Bin <- function(x){
for(i in 1:length(x))
		if(x[i]>='2.57'){
		x[i]=1
			} else{
			if(x[i]>'0'){
			x[i]=0
				} else{
				if(x[i]>='-2.57'){
				x[i]=1
					} else{
					if(x[i]>'-0'){
					x[i]=0}}
}
}
x
}


## EBDBN
## FUNCTION TO EXTRACT REGULATION INTO CYTOSCAPE FORMAT (FROM RAU)
visualize <- function (zscores, sig, type = "feedback")  {
    row <- dim(zscores)[1]
    col <- dim(zscores)[2]
    network <- matrix(0, nrow = row, ncol = col)
    network <- sign(zscores) * sig
    pos <- matrix(which(network == 1, arr.ind = TRUE)[, 2:1], ncol = 2)
    neg <- matrix(which(network == -1, arr.ind = TRUE)[, 2:1], ncol = 2)
    pos.edges <- cbind(pos, rep(1, length(which(network == 1))))
    neg.edges <- cbind(neg, rep(-1, length(which(network ==  -1))))
    edges <- rbind(pos.edges, neg.edges)
    if(type == "input") {
        edges <- cbind(paste("TF", edges[,1], sep = ""), paste("G", edges[,2], sep = ""), edges[,3])
    }
    return(edges)
}



## Bagging for Ensemble Learning
## Bootstrap aggregation (Bagging) + 1 classifier = RMSE
bagging <- function(training, testing, m=10, ite=2, meth="svmRadial",gridZ=ctl){

	Predd <- foreach(i=1:ite,.combine=cbind,.packages='caret') %dopar% {

		bagging.index <- sample(1:dim(training)[1], size=floor((nrow(training)/m)))	## vector list

		Train.me <- train(species~.,
				data=training[bagging.index,],
				method=meth,
				tuneGrid=gridZ,
				trControl=trainControl(method="repeatedcv", number=10, repeats=5
					#classProbs=T,
					#summaryFunction=twoClassSummary
					),
				#tuneLength=tune,
				preProc=c("center","scale")
				#metric="ROC"
				)

	predict(Train.me, newdata=testing)
}
	rowMeans(Predd)
}



## source : http://tinyurl.com/lo53qls
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
{
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                      mean = mean   (xx[,col], na.rm=na.rm),
                      sd   = sd     (xx[,col], na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean"=measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}



## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  require(plyr)

  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )

  # Put the subject means with original data
  data <- merge(data, data.subjMean)

  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)

  # Remove this subject mean column
  data$subjMean <- NULL

  return(data)
}




## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {

  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))

  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }

  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL

  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)

  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")

  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )

  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor

  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

}

