=============
Abbreviations
=============
FUN=function name
c()=constructor
pkname=package name
#########################################################################
#								SECTION -1-								#
#########################################################################
==================
Cheat-Sheet (2013)
==================
search()					# searching running packages
library()					# viewing the list of installed packages
installed.packages()		# printing detail of installed packages
remove.packages()			# uninstall packages
packageDescription("limma")	# one package description
install.packages("")		# install a package
require(libraryname)		#load library
detach(package:libraryname)	# unload library
source("script.R", echo=T)			# read a text file and execute its contents, echo will print the commands b4 being executed
data(datasetname, package="pkname")	# loading dataset from a package eg MASS
edit(dataset)
data(package="pkname")		# see available datasets
vignette(package="pkname")
source("http://bioconductor.org/biocLite.R")	## bioConductor
biocLite()
demo("colors")
y <- rainbow(30)			## automatic vector creator
toLatex(sessionInfo())			## extract to LaTeX the R session info
update.packages(checkBuilt=TRUE, ask=FALSE)			## update packages
system.time()				## run time of a particular script
source("http://bioconductor.org/biocLite.R")		## update packages from bioConductor
biocLite("limma")									## update packages from bioConductor
#########################################################################
#								SECTION -2-								#
#########################################################################
	## exploring an object
str()
class()
show()
names()
summary()
print() 					#if I dont want to set an object
cat() 						#like print but not on list&matrices
head(dataset)
help.start(FUN)				# viewing the supplied documentation
example(FUN)				# getting help on a function
args(FUN)
help.search(FUN)			#will output packages w/ matching FUN
help(FUN, package="nameofpackage")
RSiteSearch("FUN")			# search the web
rainbow(12)					# create 12 color palette
par(mfrow=c(3,3))			# create plot window w 3 rows and 3 cols
#########################################################################
#								SECTION -3-								#
#########################################################################	
	## listing variables and showing their structure
ls.str()
rm(list=ls())
history()
list.files()				# list files of the working directory
x <- .Last.value			# saving the result of the previous FUN (if no variable was previously defined)
traceback()					# debugging
debug(fun) ; undebug(fun)	#debuging
object.size(ls)				# ?'memory-limits'
sessionInfo()
gc()						# reallocate RAM to windows
system.time()				# benchmark a command
#########################################################################
#								SECTION -4-								#
#########################################################################		
	## Statistis
na.rm=T 					# ignore NA values when doing basic stat ie. sd(x, na.rm=T). x is a vector
seq(from=1, to=5, by=2)		# incremental sequence
seq(from=1, to=20, length.out=5)	# specify number of observations
rep(1, time=5)				# repeat value
function(x) sd(x)/mean(x)	# create a new function. see ?Control
edit(dataframe)				# edit a data frame in a cell window
na.omit(dataframe)			# remove lines w/ NAs
(dataframe, select= -col1, subset(col2 > 10))	# exclude col1 & extract lines where col2 > 10
dataframe[, -c(columns)]			# exclude columns
rbind, cbind				# combine w/ same n of col and rows
package=reshape2 & plyr		# powerful FUN for slicing & recombining dataframes
merge() 					# common col
%*%							# matmult
y %*% x						# x = ncol(y)
DF[DF$. == "x",]			# extract all rows containing x (&, |, can be added for additional options)
matrix(c(), x, y, dimnames)	# dimnames :: add names to col and rows			
dimnames(x)[[1]]			# x is a matrix, add names to 1-rows and 2-col
FUN on duplicate rows		# ddply (plyr), tapply, by, aggregate
select & subset				# select colnames ; subset is used for rows
cov.wt()					# covariance and mean calculator + they can be  weighted
paste("g",as.character(1:10),sep="")		## Add number to letter
letters[1:10]				# generate letters
gl(n,k,n*k, label=)			# n=level, k=replipaction, n*k=length of the result	(FACTOR)
#########################################################################
#								SECTION -5-								#
#########################################################################	
	## DATA
factor = is a data frame with characters and numeric
numeric = is a data frame with only numbers
character = is a data frame with only letters or " "
LETTERS[1:10]				# create a vector of letters. use in labels (plot)
{							## Load L.Rdata
setwd("C:/Dropbox/Workshop2013/Work/R/ebdbn/")
load("networks.Rdata", .GlobalEnv)
lsos(pat="")
}

#########################################################################
#								SECTION -6-								#
#########################################################################	
	## Application
linear vs non linear = linear (= the number on the axis increase equally) nonlinear (they increase nonlinearly)
residuals= the difference between the sample and the mean of the sample
fitted values = in a simple linear regression, are the point forming the linear line
R2 = variance of a sample / residual of the sample (the variance explained); if the residula is small: so landing on the fitted line, R2 is big (~1)

#########################################################################
#								SECTION -7-								#
#########################################################################
==========
Test area
==========
===========================
(function(x,y){
				z <- x^2 + y^2 ; x+y+z})(0:7,1)
===================================

bio2013 <- data.frame(c(stage),
						Growth=gm.com$Growth,
						Survival=sm.com$Mortality,
						Performance=pm$Performance)
row.names(bio2013) <- bio2013$c.stage. ; bio2013 <- bio2013[,-1]

# 3D 1.3
library(rgl)
plot3d(bio2013[,c(1,3,2)],col="red", type="n", xlab="Growth",ylab="Performance",zlab="Survival")
text3d(bio2013[,c(1,3,2)], text=rownames(bio2013))

b <- gl(4,3,12, label=c("a.E", "b.V", "c.P", "d.J"))
test <- apply(neutral, 2, function(x) {kruskal.test(x ~ b)})
dim(test)
class(test)
length(test)












