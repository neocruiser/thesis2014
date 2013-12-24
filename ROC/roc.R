#########################################################################
#								SECTION -1-								#
#########################################################################
## NETWORK ANALYSIS OF EBDBN BAYESIAN NETWORK INFERENCES

require(pROC)
roc(aSAH$outcome, aSAH$s100b)
roc(outcome ~ ndka, aSAH)

roc(outcome ~ s100b, aSAH, smooth=TRUE)

roc1 <- roc(aSAH$outcome,
aSAH$s100b, percent=TRUE,
# arguments for auc
partial.auc=c(100, 90), partial.auc.correct=TRUE,
partial.auc.focus="sens",
# arguments for ci
ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
# arguments for plot
plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
print.auc=TRUE, show.thres=TRUE)


with(aSAH, table(gender, outcome))


#########################################################################
#								SECTION -2-								#
#########################################################################


testThis <- NULL
for(i in 1:ncol(sigBinary)){
testThis <- rbind(testThis, sigBinary[sigBinary[,i] == 1,])
}

sigBinary
sigBinary[sigBinary[,1] == 0,]

rownames(testThis)
testThis2 <- testThis[! duplicated(rownames(testThis)),]

testdf <- data.frame(X1=rownames(testThis2), X2="TRUE", X3=1)
testThis_not <- sigBinary[! rownames(sigBinary) %in% rownames(testThis2),]
testdf2 <- data.frame(X1=rownames(testThis_not), X2="FALSE", X3=0)
testdf_final <- rbind(testdf, testdf2)

roc(X2 ~ X3, testdf_final)
