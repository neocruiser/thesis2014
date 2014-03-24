## updated Wednesday, 12/2013
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
visualizeNet<- function (zscores, sig, type = "feedback")  {
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
## source : (http://tinyurl.com/ljt9f78)
bagging <- function(training, testing, m=10, ite=2, meth="svmRadial",tune=5){
	Predd <- foreach(i=1:ite,.combine=cbind,.packages='caret') %dopar% {
		bagging.index <- sample(1:dim(training)[1], size=floor((nrow(training)/m)))	## vector list
		Train.me <- train(y~.,
				data=training[bagging.index,],
				method=meth,
#				tuneGrid=gridZ,
				trControl=trainControl(method="repeatedcv",number=10,repeats=5
#                                    classProbs=T,summaryFunction=twoClassSummary
                                    ),
				tuneLength=tune,
				preProc=c("center","scale")
#				metric="ROC"
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




##########################################
## Intersect and Venn Diagram Functions ##
##########################################
## Author: Thomas Girke
## Last update: March 24, 2012
## Utilities:
## (1) Venn Intersects
##     Computation of Venn intersects among 2-20 or more sample sets using the typical
##     'only in' intersect logic of Venn comparisons, such as: objects present only in
##     set A, objects present only in the intersect of A & B, etc. Due to this restrictive
##     intersect logic, the combined Venn sets contain no duplicates.
## (2) Regular Intersects
##     Computation of regular intersects among 2-20 or more sample sets using the
##     following intersect logic: objects present in the intersect of A & B, objects present
##     in the intersect of A & B & C, etc. The approach results usually in many duplications
##     of objects among the intersect sets.
## (3) Graphical Utilities
##     - Venn diagrams of 2-5 sample sets.
##     - Bar plots for the results of Venn intersect and all intersect approaches derived
##       from many samples sets.
##
## Detailed instructions for using the functions of this script are available on this page:
##     http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn
##
## Revision history:
##     March 24, 2012: fixed substring problem in plotVenn function

#######################################
## Define Generic Intersect Function ##
#######################################
## Computation of (1) Venn Intersects and (2) Regular Intersects
overLapper <- function(setlist=setlist, complexity=1:length(setlist), sep="-", cleanup=FALSE, keepdups=FALSE, type) {
	## Clean up of sample sets to minimize formatting issues
	if(cleanup==TRUE) {
		## Set all characters to upper case
		setlist <- sapply(setlist, function(x) gsub("([A-Z])", "\\U\\1", x, perl=T, ignore.case=T))
		## Remove leading and trailing spaces
		setlist <- sapply(setlist, function(x) gsub("^ {1,}| {1,}$", "", x, perl=T, ignore.case=T))
	}

	## Append object counter to retain duplicates
	if(keepdups==TRUE) {
		dupCount <- function(setlist=setlist) {
			count <- table(setlist)
			paste(rep(names(count), count), unlist(sapply(count, function(x) seq(1, x))), sep=".")
		}
		mynames <- names(setlist)
		setlist <- lapply(setlist, function(x) dupCount(x)) # lapply necessary for numeric data!
		names(setlist) <- mynames
	}

	## Create intersect matrix (removes duplicates!)
	setunion <- sort(unique(unlist(setlist)))
	setmatrix <- sapply(names(setlist), function(x) setunion %in% unique(setlist[[x]]))
	rownames(setmatrix) <- setunion
	storage.mode(setmatrix) <- "numeric"

	## Create all possible sample combinations within requested complexity levels
	labels <- names(setlist)
	allcombl <- lapply(complexity, function(x) combn(labels, m=x, simplify=FALSE))
	allcombl <- unlist(allcombl, recursive=FALSE)
	complevels <- sapply(allcombl, length)

	## Return intersect list for generated sample combinations
	if(type=="intersects") {
		OLlist <- sapply(seq(along=allcombl), function(x) setunion[rowSums(setmatrix[, rep(allcombl[[x]], 2)]) == 2 * length(allcombl[[x]])])
		names(OLlist) <- sapply(allcombl, paste, collapse=sep)
		return(list(Set_List=setlist, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Intersect_List=OLlist))
	}

	## Return Venn intersect list for generated sample combinations
	if(type=="vennsets") {
		vennSets <- function(setmatrix=setmatrix, allcombl=allcombl, index=1) {
			mycol1 <- which(colnames(setmatrix) %in% allcombl[[index]])
			mycol2 <- which(!colnames(setmatrix) %in% allcombl[[index]])
			cond1 <- rowSums(setmatrix[, rep(mycol1, 2)]) == 2 * length(mycol1)
			cond2 <- rowSums(setmatrix[, rep(mycol2, 2)]) == 0
			return(setunion[cond1 & cond2])
		}
		vennOLlist <- sapply(seq(along=allcombl), function(x) vennSets(setmatrix=setmatrix, allcombl=allcombl, index=x))
		names(vennOLlist) <- sapply(allcombl, paste, collapse=sep)
		return(list(Set_List=setlist, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Venn_List=vennOLlist))
	}
}

###########################################
## Define Venn Diagram Plotting Function ##
###########################################
vennPlot <- function(counts=counts, mymain="Venn Diagram", mysub="default", setlabels="default", yoffset=seq(0,10,by=0.34), ccol=rep(1,31), colmode=1, lcol=c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), lines=c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), mylwd=3, diacol=1, type="ellipse", ccex=1.0, lcex=1.0, sepsplit="_", ...) {
	## Enforce list structure to support multiple venn sets
	if(is.list(counts)==FALSE) {
		counts <- list(counts)
	}

	## Check for supported number of Venn counts: 3, 7, 15 and 31
	if(!length(counts[[1]]) %in%  c(3,7,15,31)) stop("Only the counts from 2-5 way venn comparisons are supported.")

        ## Function to return for a set label the index of matches in the name field of a counts object
        grepLabel <- function(label, x=names(counts[[1]])) {
                x <- strsplit(x, sepsplit)
                as.numeric(which(sapply(x, function(y) any(y==label))))
        }

	## 2-way Venn diagram
	if(length(counts[[1]])==3) {
		## Define subtitle
		if(mysub=="default") {
                        n <- names(counts[[1]])[1:2]
                        if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
			        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			        if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  }
			} else {
                                sample_counts <- sapply(n, function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
                        }
			mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), sep="")
		} else {
			mysub <- mysub
		}

		## Plot venn shapes
		symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);

		## Add counts
		for(i in seq(along=counts)) {
			olDF <- data.frame(x=c(3.1, 7.0, 5.0),
                                           y=c(6.0, 6.0, 6.0),
                                           counts=counts[[i]])
                        if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                        if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) } # For coloring several numbers per intersect differently. ccol can needs to be list to color each field differently..
		}

		## Add sample labels
		if(length(setlabels)==1 & setlabels[1]=="default") {
			setlabels <- names(counts[[1]][1:2])
		} else {
			setlabels <- setlabels
		}
		text(c(2.0, 8.0), c(8.8, 8.8), labels=setlabels, col=lcol, cex=lcex, ...)
	}

	## 3-way Venn diagram
	if(length(counts[[1]])==7) {
		## Define subtitle
		if(mysub=="default") {
                        n <- names(counts[[1]])[1:3]
                        if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
			        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			        if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  }
			} else {
			        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			}
                        mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), sep="")
		} else {
			mysub <- mysub
		}

		## Plot venn shapes
		symbols(x=c(4, 6, 5), y=c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=FALSE, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", fg=lines, ...)

		## Add counts
		for(i in seq(along=counts)) {
			olDF <- data.frame(x=c(3.0, 7.0, 5.0, 5.0, 3.8, 6.3, 5.0),
                                           y=c(6.5, 6.5, 3.0, 7.0, 4.6, 4.6, 5.3),
                                           counts=counts[[i]])
	        	 if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                         if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }

		}

                ## Add sample labels
		if(length(setlabels)==1 & setlabels[1]=="default") {
			setlabels <- names(counts[[1]][1:3])
		} else {
			setlabels <- setlabels
		}
		text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), labels=setlabels, col=lcol, cex=lcex, ...)
	}

	## 4-way Venn diagram with ellipses
	if(length(counts[[1]])==15 & type=="ellipse") {
		## Define subtitle
		if(mysub=="default") {
                        n <- names(counts[[1]])[1:4]
                        if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
			        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			        if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  }
			} else {
			        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			}
                        mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
		} else {
			mysub <- mysub
		}

		## Plot ellipse
		plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
			angles <- (0:segments) * 2 * pi/segments
			rotate <- rotate*pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
			ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])
			plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
		}
		## Plot ellipse as 4-way venn diagram
		ellipseVenn <- function(...) {
			split.screen(c(1,1))
			plotellipse(center=c(3.5,3.6), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=mymain, sub=mysub, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(4.7,4.4), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(5.3,4.4), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.5,3.6), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, lwd=mylwd, ...)
			## Add counts
			for(i in seq(along=counts)) {
				olDF <- data.frame(x=c(1.5, 3.5, 6.5, 8.5, 2.9, 3.1, 5.0, 5.0, 6.9, 7.1, 3.6, 5.8, 4.2, 6.4, 5.0),
                                                   y=c(4.8, 7.2, 7.2, 4.8, 5.9, 2.2, 0.7, 6.0, 2.2, 5.9, 4.0, 1.4, 1.4, 4.0, 2.8),
                                                   counts=counts[[i]])
	        	        if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                                if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }
			}
			## Add sample labels
			if(length(setlabels)==1 & setlabels[1]=="default") {
				setlabels <- names(counts[[1]][1:4])
			} else {
				setlabels <- setlabels
			}
			text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 8.3, 8.3, 7.3), labels=setlabels, col=lcol, cex=lcex, ...)
			close.screen(all=TRUE)
		}
		ellipseVenn(...)
	}

	## 4-way Venn diagram with circles (pseudo-venn diagram that misses two overlap sectors)
	if(length(counts[[1]])==15 & type=="circle") {
		## Define subtitle
		if(mysub=="default") {
                        n <- names(counts[[1]])[1:4]
                        if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
			        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			        if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  }
			} else {
			        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			}
                        mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
		} else {
			mysub <- mysub
		}

		## Plot venn shapes
		symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=FALSE, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", fg=lines, ...)

		## Add counts
		for(i in seq(along=counts)) {
		        olDF <- data.frame(x=c(3.0, 6.5, 3.0, 6.5, 4.8, 3.0, 4.8, 4.8, 6.5, 4.8, 3.9, 5.7, 3.9, 5.7, 4.8),
                                           y=c(7.2, 7.2, 3.2, 3.2, 7.2, 5.2, 0.4, 0.4, 5.2, 3.2, 6.3, 6.3, 4.2, 4.2, 5.2),
                                           counts=counts[[i]])
			 if(colmode==1) { text(olDF$x[-c(7,8)], olDF$y[-c(7,8)] + yoffset[i], olDF$counts[-c(7,8)], col=ccol, cex=ccex, ...) } # rows 14-15 of olDF are printed in next step
			 if(colmode==2) { text(olDF$x[-c(7,8)], olDF$y[-c(7,8)] + yoffset[i], olDF$counts[-c(7,8)], col=ccol[[i]], cex=ccex[i], ...) }
			 text(c(4.8), c(0.8) + yoffset[i], paste("Only in ", names(counts[[1]][1]), " & ", names(counts[[1]][4]), ": ", olDF$counts[7], "; Only in ", names(counts[[1]][2]), " & ", names(counts[[1]][3]), ": ", olDF$counts[8], sep=""), col=diacol, cex=ccex, ...)
		}

                ## Add sample labels
			if(length(setlabels)==1 & setlabels[1]=="default") {
				setlabels <- names(counts[[1]][1:4])
			} else {
				setlabels <- setlabels
			}
		text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), labels=setlabels, col=lcol, cex=lcex, ...)
	}

	## 5-way Venn diagram
	if(length(counts[[1]])==31) {
		## Define subtitle
		if(mysub=="default") {
                        n <- names(counts[[1]])[1:5]
                        if(!all(rowSums(sapply(n, function(x) sapply(n, function(y) grepl(y, x)))) == 1)) { # Checks if one or more set labels are substrings of one another
			        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			        if(sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, names(counts[[1]][-c(1:length(n))])))) { sample_counts <- rep("?", length(n)); warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")  }
			} else {
			        sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, names(counts[[1]]))]))
			}
                        mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), paste("; S5 =", sample_counts[5]), sep="")
		} else {
			mysub <- mysub
		}

		## Plot ellipse
		plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
			angles <- (0:segments) * 2 * pi/segments
			rotate <- rotate*pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
			ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])
			plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
		}
		## Plot ellipse as 5-way venn diagram
		ellipseVenn <- function(...) {
			split.screen(c(1,1))
			screen(1, new=FALSE)
			plotellipse(center=c(4.83,6.2), radius=c(1.43,4.11), rotate=0, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=mymain, sub=mysub, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.25,5.4), radius=c(1.7,3.6), rotate=66, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.1,3.5), radius=c(1.55,3.9), rotate=150, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(4.48,3.15), radius=c(1.55,3.92), rotate=210, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(3.7,4.8), radius=c(1.7,3.6), rotate=293.5, segments=360, xlab="", ylab="", col=lines[5], axes=FALSE, lwd=mylwd, ...)

			## Add counts
			for(i in seq(along=counts)) {
				olDF <- data.frame(x=c(4.85, 8.0, 7.1, 3.5, 2.0, 5.90, 4.4, 4.60, 3.60, 7.1, 6.5, 3.2, 5.4, 6.65, 3.40, 5.00, 6.02, 3.60, 5.20, 4.03, 4.20, 6.45, 6.8, 3.39, 6.03, 5.74, 4.15, 3.95, 5.2, 6.40, 5.1),
                                                   y=c(8.30, 6.2, 1.9, 1.6, 5.4, 6.85, 6.6, 2.45, 6.40, 4.3, 6.0, 4.6, 2.1, 3.40, 3.25, 6.43, 6.38, 5.10, 2.49, 6.25, 3.08, 5.30, 4.0, 3.80, 3.20, 5.95, 5.75, 3.75, 3.0, 4.50, 4.6),
					counts=counts[[i]])
	        	        if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                                if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }
			}
			## Add sample labels
			if(length(setlabels)==1 & setlabels[1]=="default") {
				setlabels <- names(counts[[1]][1:5])
			} else {
				setlabels <- setlabels
			}
			text(c(5.7, 7.9, 8.5, 4.2, 0.8), c(9.9, 7.9, 1.9, 0.0, 7.3), adj=c(0, 0.5), labels=setlabels, col=lcol, cex=lcex, ...)
			close.screen(all=TRUE)
		}
		ellipseVenn(...)
	}
}

##############################
## Define Bar Plot Function ##
##############################
## Plots the counts of Venn/regular intersects generated by the overLapper function
olBarplot <- function(OLlist=OLlist, mycol="default", margins=c(6, 10, 3, 2), mincount=0, mysub="default", ...) {
	## Generate counts and allow lower limit
	counts <- sapply(OLlist[[4]], length)
	mylogical <- counts >= mincount
	counts <- counts[mylogical]

	## Color bars by default by complexity levels
	if(mycol=="default") {
		mycol <- OLlist$Complexity_Levels
		mycol <- mycol[mylogical]
	} else {
		mycol <- mycol
	}

	## Define subtitle
	if(mysub=="default") {
		mysub <- paste("Min Count:", mincount)
	} else {
		mysub <- mysub
	}

	## Generate bar plot with defined margins
	par(mar=margins) # Define margins to allow long labels
	barplot(counts, col=mycol, sub=mysub, ...)
	par(mar=c(5, 4, 4, 2) + 0.1) # Set margins back to default
}

regfit <- function(y,x,...){
print(summary(...))
fit <- lm(y~x)
par(mfrow = c(2,2))
plot(fit)
print(summary(fit))
}
regfit.int <- function(y,x,z,...){
print(summary(...))
fit <- lm(y~x*z)
par(mfrow = c(2,2))
plot(fit)
print(summary(fit))
}
## Linear regression wrap up functions

predict.regsubsets <- function(object, newdata, id, ...){
form <- as.formula(object$call[[2]])
mat <- model.matrix(form, newdata)
coefi <- coef(object, id=id)
mat[,names(coefi)]%*%coefi
}
## source 06-3-R for forward and backward subset selection form Stat learning ch6


rocplot <- function(pred, truth, ...){
predob <- prediction(pred, truth)
perf <- performance(predob, "tpr", "fpr")
plot(perf,...)}
## plot ROC curves (p365 ch9 Stat leanring)

Cols <- function(vec){
cols <- rainbow(length(unique(vec)))
return(cols[as.numeric(as.factor(vec))])
}
## return colours (used in PCA ann2.R)

set.var <- function(dat,n,b){
locus.var <- apply(t(dat), n, var)
hi.var <- order(abs(locus.var), decreasing = T)[1:b]
cat("\nRange of the variance:",range(locus.var),"\n")
cat("\nNumber of selected high-variance genes:",length(hi.var),"\n")
return(hi.x <- dat[,hi.var])
}
## Unsupervised gene selection based on high variance

set.cor <- function(x,cutoff){
require(caret)
cor.x <- cor(hi.x)
print(summary(cor.x[upper.tri(cor.x)]))
cor.hi.x <- findCorrelation(cor.x, cutoff=cutoff)
return(hi.x <- hi.x[, -cor.hi.x])
}
## Unsupervised gene selection based on low correlation

range_correlation <- function(dataY,n=1,t=1,data,method){
resultats <- solutions(dataY)[[n]][t,]
test=NULL
for(i in 1:length(resultats)){test <- c(test,correlate(data[,n],data[,resultats[i]],method=method)$statistic)}
cat("\n","The range of the correlations:")
return(range(test))
}
## require(mRMRe). Compute correlation between selected features (feature selection Paper3) then the range

locusRMR <- function(dataY,dataX,n,t){
print(apply(solutions(dataY)[[n]], 2, function(x,y) { return(y[x]) }, y=featureNames(dataX))[t,])
}
## require(mRMRe). returns locus names


model.reg <- function(dat,train,test,method,folds=10,rep=5,tune){
    trainCtrl <- trainControl(method="repeatedcv",number=folds, repeats=rep)	## Regression
    lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                              method=method,
                                              trControl= trainCtrl,
                                              preProc=c("center","scale"),
                                              tuneLength=tune ))
    ploted <- plot(modelTrain)
    Predd <- predict(modelTrain, newdata=dat[test,], type="raw")
    rmse <- mean((Predd - y[test])^2)		## Test set MSE for regression
    output <- list(ploted,TimeLapsed=lapsed, Hyperparameters=modelTrain$bestTune, RMSE=rmse)
    return(output)
}
## require(caret). Regression. Train model and test on independant test set return(RMSE)

modelTune.reg <- function(dat,train,test,method,folds=10,rep=5,tune,ctl){
    trainCtrl <- trainControl(method="repeatedcv",number=folds, repeats=rep)	## Regression
    lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                              method=method,
                                              trControl= trainCtrl,
                                              preProc=c("center","scale"),
                                              tuneGrid=ctl,
                                              tuneLength=tune ))
    ploted <- plot(modelTrain)
    Predd <- predict(modelTrain, newdata=dat[test,], type="raw")
    rmse <- mean((Predd - y[test])^2)		## Test set MSE for regression
    output <- list(ploted,TimeLapsed=lapsed, Hyperparameters=modelTrain$bestTune, RMSE=rmse)
    return(output)
}
## require(caret). Regression. Uses GRID for HYPERPARAMETER tuning. Train model and test on independant test set return(RMSE)

model.clas <- function(dat,train,test,method,folds=10,rep=5){
    trainCtrl <- trainControl(method="repeatedcv",number=folds, repeats=rep, , classProbs=T,summaryFunction=defaultSummary) ## classification
    lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                              method=method,
                                              trControl= trainCtrl,
                                              preProc=c("center","scale"),
                                              tuneLength=tune,
                                              metric="ROC"))
    ploted <- plot(modelTrain)
    Predd <- predict(modelTrain, newdata=dat[test,], type="raw")
    conf.m <- confusionMatrix(data=Predd, dat[test,1])	## confusion matrix for classification
    Probs <- predict(modelTrain, newdata=dat[test,], type="prob")
    output <- list(ploted,TimeLapsed=lapsed, Hyperparameters=modelTrain$bestTune, ConfusionMatrix=conf.m,Probabilities=Probs)
    return(output)
}
## require(caret). Classification. Train model and test on independant test set return(classification error)

modelTune.clas <- function(dat,train,test,method,folds=10,rep=5,tune){
    trainCtrl <- trainControl(method="repeatedcv",number=folds, repeats=rep, , classProbs=T,summaryFunction=defaultSummary) ## classification
    lapsed <- system.time(modelTrain <- train(y~., data=dat[train,],
                                              method=method,
                                              trControl= trainCtrl,
                                              preProc=c("center","scale"),
                                              tuneGrid=ctl,
                                              tuneLength=tune,
                                              metric="ROC"))
    ploted <- plot(modelTrain)
    Predd <- predict(modelTrain, newdata=dat[test,], type="raw")
    conf.m <- confusionMatrix(data=Predd, dat[test,1])	## confusion matrix for classification
    Probs <- predict(modelTrain, newdata=dat[test,], type="prob")
    output <- list(ploted,TimeLapsed=lapsed, Hyperparameters=modelTrain$bestTune, ConfusionMatrix=conf.m,Probabilities=Probs)
    return(output)
}
## require(caret). Classification. GRID search HYPERPARAMETERS tuning. Train model and test on independant test set return(classification error)









