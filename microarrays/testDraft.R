
names(MA)

df2 <- MA2[MA2$SystematicName == "gle1", 3:5]

results1.4 <- subset(results, results$oeuf2 > 0.2)
df8 <- MA[MA$M > 1, ]
results[1:10,]

ndf <- results[results$X2>1, ]

class(MA)
str(MA)
#results contains Factors. i should convert to Vectors to SORT w/ subset

str(results)
class(results)
unclass(results)
results.c <- as.character(results)
names(results.c)
class(results.n)
str(results.c)
ndf <- subset(results.n, results.ve$oeuf2 > 0.2)
class(results.v)

oeuf1.r <- results[,2]
oeuf1.r
class(oeuf1.r)
str(oeuf1.r)
unclass(oeuf1.r)
as.numeric(oeuf1.r)
as.character(oeuf1.r)
oeuf1.rdf <- data.frame(oeuf1.r)
names(oeuf1.rdf) <- c("oeuf1")
names(oeuf1.rdf)
results0.2 <- subset(oeuf1.rdf, oeuf1.rdf$oeuf1 > 0.2)
rm(oeuf1.r)