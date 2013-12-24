### Supplementary material with the paper
### Interpretation of ANOVA models for microarray data using PCA
### J.R. de Haan, R. Wehrens, S. Bauerschmidt, E. Piek,
### R.C. van Schaik and L.M.C. Buydens
### Bioinformatics (2006)
###
### Please cite this paper when you use this code in a publication.
###
### Written by J.R. de Haan, December 18, 2006
###

### Step1: a synthetic dataset of 500 genes is generated with 5 classes
###	1 unresponsive genes (300 genes)
###	2 constant genes (50 genes)
###	3 profile 1 (50 genes)
###	4 profile 2 (50 genes)
###	5 profile 3 (50 genes)

#generate synthetic dataset with similar dimensions:
# 500 genes, 3 replicates, 10 timepoints, 4 treatments
X <- array(0, c(500, 3, 10, 4))
labs.synth <- c(rep(1, 300), rep(2, 50), rep(3, 50), rep(4, 50), rep(5, 50))
gnames <- cbind(labs.synth, labs.synth)
#print(dim(gnames))
gnames[1:300,2] <- "A"
gnames[301:350,2] <- "B"
gnames[351:400,2] <- "C"
gnames[401:450,2] <- "D"
gnames[451:500,2] <- "E"

### generate 300 "noise" genes with expressions slightly larger than
### the detection limit (class 1)
X[labs.synth==1,1,,] <- rnorm(length(X[labs.synth==1,1,,]), mean=50, sd=40)
X[labs.synth==1,2,,] <- X[labs.synth==1,1,,] +
  rnorm(length(X[labs.synth==1,1,,]), mean=0, sd=10)
X[labs.synth==1,3,,] <- X[labs.synth==1,1,,] +
  rnorm(length(X[labs.synth==1,1,,]), mean=0, sd=10)

# generate 50 stable genes at two levels (class 2)
X[301:325,1,,] <- rnorm(length(X[301:325,1,,]), mean=1500, sd=40)
X[301:325,2,,] <- X[301:325,1,,] + rnorm(length(X[301:325,1,,]), mean=0, sd=10)
X[301:325,3,,] <- X[301:325,1,,] + rnorm(length(X[301:325,1,,]), mean=0, sd=10)

X[326:350,1,,] <- rnorm(length(X[326:350,1,,]), mean=3000, sd=40)
X[326:350,2,,] <- X[326:350,1,,] + rnorm(length(X[326:350,1,,]), mean=0, sd=10)
X[326:350,3,,] <- X[326:350,1,,] + rnorm(length(X[326:350,1,,]), mean=0, sd=10)

# generate50 genes with profile 1 (class 3)
increase.range <- matrix(rep(1:50, 10), ncol=10, byrow=FALSE)
profA3 <- matrix(rep(c(10, 60, 110, 150, 150, 150, 150, 150, 150, 150) , 50),
                 ncol=10, byrow=TRUE) * increase.range
X[351:400,1,,1] <- profA3 + rnorm(length(profA3), mean=0, sd=40)
profB3 <- matrix(rep(c(10, 100, 220, 280, 280, 280, 280, 280, 280, 280), 50),
                 ncol=10, byrow=TRUE) * increase.range
X[351:400,1,1:10,2] <- profB3 + rnorm(length(profA3), mean=0, sd=40)
profC3 <- matrix(rep(c(10, 120, 300, 300, 280, 280, 280, 280, 280, 280), 50),
                 ncol=10, byrow=TRUE) * increase.range
X[351:400,1,1:10,3] <- profC3 + rnorm(length(profA3), mean=0, sd=40)
profD3 <- matrix(rep(c(100, 75, 50, 50, 50, 50, 50, 50, 75, 100), 50), ncol=10,
                 byrow=TRUE)
X[351:400,1,1:10,4] <- profD3 + rnorm(length(profA3), mean=0, sd=40)
#again replicates
X[351:400,2,,] <- X[351:400,1,,] + rnorm(length(X[351:400,2,,]), mean=0, sd=10)
X[351:400,3,,] <- X[351:400,1,,] + rnorm(length(X[351:400,3,,]), mean=0, sd=10)

# generate50 genes with profile 2 (class 4)
increase.range <- matrix(rep(1:50, 10), ncol=10, byrow=FALSE)
profA4 <- matrix(rep(c(10, 60, 110, 150, 125, 100, 75, 50, 50, 50) , 50),
                 ncol=10, byrow=TRUE) * increase.range
X[401:450,1,,1] <- profA4 + rnorm(length(profA4), mean=0, sd=40)
profB4 <- matrix(rep(c(10, 100, 220, 280, 200, 150, 100, 50, 50, 50), 50),
                 ncol=10, byrow=TRUE) * increase.range
X[401:450,1,1:10,2] <- profB4 + rnorm(length(profA4), mean=0, sd=40)
profC4 <- matrix(rep(c(10, 150, 300, 220, 150, 100, 50, 50, 50, 50), 50),
                 ncol=10, byrow=TRUE) * increase.range
X[401:450,1,1:10,3] <- profC4 + rnorm(length(profA4), mean=0, sd=40)
profD4 <- matrix(rep(c(150, 100, 50, 50, 75, 75, 75, 100, 100, 100), 50),
                 ncol=10, byrow=TRUE)
X[401:450,1,1:10,4] <- profD4 + rnorm(length(profA4), mean=0, sd=40)
#again replicates
X[401:450,2,,] <- X[401:450,1,,] + rnorm(length(X[401:450,2,,]), mean=0, sd=10)
X[401:450,3,,] <- X[401:450,1,,] + rnorm(length(X[401:450,3,,]), mean=0, sd=10)

# generate50 genes with profile 3 (class 5)
increase.range <- matrix(rep(1:25, 20), ncol=10, byrow=FALSE)
profA4 <- matrix(rep((200 - c(10, 60, 110, 150, 125, 100, 75, 50, 50, 50)), 50),
                 ncol=10, byrow=TRUE) * increase.range
X[451:500,1,,1] <- profA4 + rnorm(length(profA4), mean=0, sd=40)
profB4 <- matrix(rep((200 - c(10, 100, 180, 200, 200, 150, 100, 50, 50, 50)), 50),
                 ncol=10, byrow=TRUE) * increase.range
X[451:500,1,1:10,2] <- profB4 + rnorm(length(profA4), mean=0, sd=40)
profC4 <- matrix(rep((200 - c(10, 150, 200, 180, 150, 100, 50, 50, 50, 50)), 50),
                 ncol=10, byrow=TRUE) * increase.range
X[451:500,1,1:10,3] <- profC4 + rnorm(length(profA4), mean=0, sd=40)
profD4 <- matrix(rep((200 - c(150, 100, 50, 50, 75, 75, 75, 100, 100, 100)), 50),
                 ncol=10, byrow=TRUE)
X[451:500,1,1:10,4] <- profD4 + rnorm(length(profA4), mean=0, sd=40)
#again replicates
X[451:500,2,,] <- X[451:500,1,,] + rnorm(length(X[451:500,2,,]), mean=0, sd=10)
X[451:500,3,,] <- X[451:500,1,,] + rnorm(length(X[451:500,3,,]), mean=0, sd=10)


# Step 2: Now the effects for different factors in the ANOVA model
# can be calculated:

# subtraction of the general mean
x <- X - mean(X, na.rm=TRUE)
tpoints <- c(1, 3, 6, 12,c( 24*c(1, 2, 3, 5, 8, 12)))
nrgenes <- dim(x)[1]

# calculation of the three main effects
cat("calculating main effects\n")
timemeans <- apply(x, 3, mean, na.rm=TRUE)
treatmeans <- apply(x, 4, mean, na.rm=TRUE)
genemeans <- apply(x, 1, mean, na.rm=TRUE)

#par(mfrow=c(2, 3))


# calculation of the interaction effects
# interaction time-treatment
cat("calculating interaction time-treatment\n")
mean.ti.tr <- apply(x, c(3,4), mean, na.rm=TRUE)
#print(dim(mean.ti.tr))
time1.m <- matrix(rep(timemeans, 4), nrow=10, byrow=FALSE)
tr1.m <- matrix(rep(treatmeans, 10), nrow=10, byrow=TRUE)
int.ti.tr <- (mean.ti.tr - time1.m) - tr1.m

# interaction time-gene
cat("calculating interaction time-gene\n")
mean.ti.gene <- apply(x, c(1,3), mean, na.rm=TRUE)
#print(dim(mean.ti.gene))
time2.m <- matrix(rep(timemeans, dim(x)[1]), nrow=dim(x)[1], byrow=TRUE)
gene2.m <- matrix(rep(genemeans, 10), nrow=dim(x)[1], byrow=FALSE)
int.ti.gene <- (mean.ti.gene - time2.m) - gene2.m

# interaction gene-treatment
cat("calculating interaction gene-treatment\n")
mean.gene.tr <- apply(x, c(1,4), mean, na.rm=TRUE)
#print(dim(mean.gene.tr))
gene3.m <- matrix(rep(genemeans, 4), ncol=4, byrow=FALSE)
tr3.m <- matrix(rep(treatmeans, dim(x)[1]), ncol=4, byrow=TRUE)
int.gene.tr <- (mean.gene.tr - gene3.m) - tr3.m

# calculation of 3 factor interaction
cat("calculating 3 factor interaction\n")
mean.gene.time.tr <- apply(x, c(1, 3, 4), mean, na.rm=TRUE)
#print(dim(mean.gene.time.tr))

ar.ti.tr <- array(0, c(nrgenes, 10, 4))
for (i in 1:nrgenes){ar.ti.tr[i,,] <- int.ti.tr}

ar.ti.gene <- array(0, c(nrgenes, 10, 4))
for (i in 1:4){ar.ti.gene[,,i] <- int.ti.gene}

ar.gene.tr <- array(0, c(nrgenes, 10, 4))
for (i in 1:10){ar.gene.tr[,i,] <- int.gene.tr}

ar.ti <- array(0, c(nrgenes, 10, 4))
for (i in 1:4){ar.ti[,,i] <- time2.m}

ar.tr <- array(0, c(nrgenes, 10, 4))
for (i in 1:10){ar.tr[,i,] <- tr3.m}

ar.gene <- array(0, c(nrgenes, 10, 4))
for (i in 1:4){ar.gene[,,i] <- gene2.m}

int.gene.time.tr <- mean.gene.time.tr - ar.ti - ar.tr - ar.gene -
  ar.ti.tr - ar.ti.gene - ar.gene.tr
imat.gtt <- cbind(int.gene.time.tr[,,1],
                  int.gene.time.tr[,,2],
                  int.gene.time.tr[,,3],
                  int.gene.time.tr[,,4])

### calculation of error term
error.term1 <- abs(sweep(x, c(1, 3, 4), mean.gene.time.tr))
cell.error <- apply(error.term1, c(1, 3, 4), mean, na.rm=TRUE)
mncn.error <- scale(cbind(cell.error[,,1], cell.error[,,2],
                          cell.error[,,3], cell.error[,,4]), scale=FALSE)

### Step 3: The results of the model can now be inspected with
### different plots (PCA for the interactions)
#plot(timemeans, type="l")
#plot(treatmeans, type="l")
#plot(genemeans, type="l")

#source("biplot.R")
#jbiplot(1, 2, int.ti.gene, gnames[,2], as.character(tpoints), rep(2, 10))
#jbiplot(1, 2, int.gene.tr, gnames[,2], c("TR1", "TR2", "TR3", "UNT"), 2:5)
#jbiplot(1, 2, int.ti.tr, as.character(tpoints), c("TR1", "TR2", "TR3", "UNT"), rep(2, 4))
#jbiplot(1, 2, imat.gtt, gnames[,2], as.character(rep(tpoints, 4)), sort(rep(2:5, 10)))
