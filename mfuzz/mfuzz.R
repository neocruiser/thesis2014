

#########################################################################
#								SECTION -1-								#
#########################################################################
## Functional enrichment analysis
## BLAST2GO


## Extract locus from blast2GO
## Extract logs
## Convert to expressionSet
## Run soft clustering






## create a minimal ExpressionSet & run MFUZZ
## ********************* GO Term *****************
## SET i and j
i=15 ; j=1.5
{
library(Biobase)
minimalSet_ns <- ExpressionSet(assayData=average)
## Soft clustering and visualisation
library(Mfuzz)
## Standardization (relative expression vs absolute ampliture)
minimalSet <- standardise(minimalSet_ns)
## Fuzzifier m testing
cat("\n\n 'm' estimate for genes",mestimate(minimalSet),"\n")
## RUN soft clustering
cl <- mfuzz(minimalSet,c=i,m=j)
mfuzz.plot(minimalSet,cl=cl,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
}
## coupling between clusters (high values between clusters = smilar pattern)
{
O <- overlap(cl)
over_plot <- overlap.plot(cl, over=O, thres=0.05)
cat("\n","\n","If coupling is large, cluster patterns are similar","\n")
over_plot
## Analysis of global clustering structures
cl2 <- mfuzz(minimalSet,c=i,m=j)
mfuzz.plot(minimalSet,cl=cl2,mfrow=c(3,3),time.labels=c("E","T","V","P","J"))
O2 <- overlap(cl2)
over_plot <- overlap.plot(cl2, over=O2,P=over_plot, thres=0.05)
}
## ADDITIONAL OPTIONS (EXTRACT CLUSTERS)
## See available components & extract GENES
## **************Subset C2
{
cl
## Plot the center of the third cluser
********* History > Recording (after flashing par())
par(mfrow=c(2,i/2))
for(z in 1:i) plot(cl[[1]][z,],type="l",ylab="Expression")
## Extract the first 10 membership of the first cluster
cl[[4]][1:10,1]
## Extract genes forming the alpha cores of soft clusters
acore_list <- acore(minimalSet, cl=cl, min.acore=0.7)
str(acore_list)
for(i in 1:length(acore_list)) print(nrow(acore_list[[i]]))
}