#=============================================================================
#=============================================================================
#Full compilation without images, plots and volcanoes
rm(MA2, MA, RG)

## START

library (limma)
library (statmod)
setwd("C:/Users/Sleiman Bassim/Dropbox/Workshop 2013/Work environment/Agilent/5- Workplace/")
targets <- readTargets("targetsO.txt")
RG <- read.maimages(targets, source="agilent")
RG <- backgroundCorrect(RG, method="normexp", offset=20, normexp.method="mle")
#===============================================================================
#test (not working) p92 limma userguide
#process to combine 2 of the same probres into one, either a correlation or an average
corfit <- duplicateCorrelation(MA, ndups= 2, spacing="rows")
corfit$consensus.correlation

#not fully working, maybe bcause RG doesnt contain spot quality weights
#nonDE contains only TERA probes, non differentially expressed
#voir smyth 2003
layout <- list(ngrid.r=1, ngrid.c=1, nspot.r=384, nspot.c=164)
MA <- normalizeWithinArrays(RG, layout)
nonDE <- grep("TERA", RG$genes[,"SystematicName"], ignore.case=FALSE, value=TRUE, invert=FALSE)
MA <- normalizeWithinArrays(RG, layout, method="composite", controlspots=nonDE)

#working
#upweighting control spots of non deferential expressions between arrays, using ControlType column
w <- modifyWeights(array(1,dim(RG)), 
                   RG$genes$ControlType, 
                   c("0","1","-1"), 
                   c(2,0,0))


#remove bacterial spots
#notepad++ search ^[^>]*$ and replace with 'nothing' (removes all sequences, keeps the heading of every fasta sequence)
#notepad++ search \n\n (extended mode) and replace w/ \n (removes all blanck lines)
#notepad++ search ^[>].*?(?=L) and replace w/ 'nothing' (removes every string before 'Locus/Node')
#notepad++ search \|.*$ and replace w/ 'nothing (removes every string after 'Locus/Node')
#R remove duplicate entries (rows), keep unique in a 'dbs' data frame
##RUN ONCE (just to create CUSTbacterialSpots file w/ 6984 spots)
#START
#remove bacterial spots from librairy
bacterialspots <- readTargets("bacterialSpots.txt")
#which row is a repeat?
duplicated(bacterialspots)
#show repeat entries
bacterialspots[duplicated(bacterialspots),]
#show unique repeat entries (w/o counting the repeats)
unique(bacterialspots[duplicated(bacterialspots),])
#Original bacterial data (from MG-RAST) w/ Locus repeats removed
dbs <- bacterialspots[!duplicated(bacterialspots),]   #bacterial locus cleaned
gal <- readGAL()
cust <- gal$ID
locus <- gal$TopHit
id <- cbind(cust, locus)
#create a text file w/ all spots under 2 columns (CUST and Locus)
id <- data.frame(id)
write.table(id, "id.txt", sep="\t")     #all locus raw
#remove quotation marks (") from id.txt notepad++
#remove \|.*$ post-Locus string form id.txt notepad++ 
#create a data frame contaning only CUSTs of bacterialspots Locus's
id.t <- readTargets("id.txt")       #all locus cleaned
names(id.t)
#remove nonbacterial spots from data frame using bacterial spots (dbs)
#note: I should convert id.txt into a data frame, otherwise I can't remove spots
# %in% (search for in ?match)
id.dbs <- id.t[id.t$locus %in% c(dbs), ]    #CUST and Locus bacteria raw
write.table(id.dbs, "CUSTbacterialSpots.txt", sep="\t")    #CUST bacteria cleaned
#remove everything except CUSTs notepad++
#notepad++ search ^.*?(?=C) and replace w/ 'nothing' (removes every string before 'CUST')
#notepad++ search \t.*$ and replace w/ 'nothing' (removes every string after 'CUST')
CUST <- read.table("CUSTbacterialSpots.txt", sep="\t")    #create data frame
#use the input in CUSTbacterialspots.txt as line code in R (old fashioned)
#add quotation marks
#remove end of line w/ notepad++ search: [\r\n]+ and replace: 'a comma and space'
#add the end of R code line w/ notepad++ search: 4", $ or 4",$ and replace: / 4"\), \]
#add the begining of the R code line w/ notepad++ search: ^"C and replace: MA <- MA[! MA$genes$ProbeName %in% c\("C

duplicated(CUST)                   #remove duplicates
CUST.t <- CUST[!duplicated(CUST),]
CUST.tt <- CUST.t[!CUST.t=="cust"]    #remove first line 'cust'
#not working
x <- read.table("CUST.tt.txt")
CUST.tt <- data.frame(x)
MA.q <- MA[!MA$genes$ProbeName %in% c(CUST.tt$x),]
dim(MA.q)
rm(MA.q)
##END

#test
gal <- readGAL()
layout <- getLayout(gal, guessdups=TRUE)
layout
getSpacing("topbottom",list(ngrid.r=1,ngrid.c=1,nspot.r=384,nspot.c=164))
printorder(list(ngrid.r=1,ngrid.c=1,nspot.r=384,nspot.c=164))

#===============================================================================
MA <- normalizeWithinArrays(RG, method="printtiploess")
MA <- normalizeBetweenArrays(MA, method="scale")
MA <- MA[MA$genes$ControlType=="0" & !MA$genes$SystematicName=="---NA---" & !MA$genes$SystematicName=="titration", ]

#add missing columns and create gene table
results1 <- MA$genes$ProbeName
results2 <- MA$M
results3 <- cbind(results1, results2)
colnames(results3) <- MA$genes$names
results <- data.frame(results3)
write.table(results, "results2.txt", sep="\t")

## END

rm(list=ls())
#=================================================================================================
#=================================================================================================
#phase1
#launch limma library, then change wrking directory
library (limma)
setwd("C:/Users/Sleiman Bassim/Dropbox/Workshop 2013/Work environment/Agilent/5- Workplace/")
getwd()

#phase2
#create object (targets)
targets <- readTargets("targetsO.txt")
RG <- read.maimages(targets, source="agilent")
plotDensities(RG)

names(RG$genes)
#(trial but not ncessary)reading probe annotations from a GAL file (source=agilent)
#refer to smyth 2003 p7
RG$genes <- readGAL()
names(RG$genes)
RG$printer <- getLayout(RG$genes)
names(RG$printer)
show(RG$printer)

#phase3
#background correction (previously used "auto" w/ offset =1)
#see smyth 2003
RG <- backgroundCorrect(RG, method="normexp", offset=15, normexp.method="mle")
plotDensities(RG)
plotMA(RG)

#test (not working)
#cant use it because Agilent Arrays don't have print-tip
#view p28 limma user guide
#layout <- list(ngrid.r=328, ngrid.c=192, nspot.r=8, nspot.c=8)
#MA <- normalizeWithinArrays(RG)

#phase4
#normalize within & between arrays. MA holds the normalized data
#within array normalization (2 channel necessary)
#global loess doesn't take subarray variation (avoid)
# view Oshlack 2007 for citation
MA <- normalizeWithinArrays(RG, method="loess")
plotDensities(MA)
boxplot(MA$M~col(MA$M),names=colnames(MA$M))
plotMA(MA)

#phase4'
#remove control probres
MA2 <- MA[MA$genes$ControlType=="0", ]
#check if Min. & Max. of controlType are 0, not 1 & -1
summary(MA2$genes)
#check if ProbeName list starts w/ experimental genes, not control probes
print(MA2)

#phase4''
#remove TERA probes from MA2 (manually)
MA3 <- MA2[!MA2$genes$SystematicName=="TERA00955", ]
MA4 <- MA3[!MA3$genes$SystematicName=="TERA01252", ]
MA5 <- MA4[!MA4$genes$SystematicName=="TERA01866", ]
MA6 <- MA5[!MA5$genes$SystematicName=="TERA01890", ]
MA7 <- MA6[!MA6$genes$SystematicName=="TERA02153", ]
MA8 <- MA7[!MA7$genes$SystematicName=="TERA02352", ]
MA9 <- MA8[!MA8$genes$SystematicName=="TERA02367", ]
MA10 <- MA9[!MA9$genes$SystematicName=="TERA02369", ]
MA11 <- MA10[!MA10$genes$SystematicName=="TERA02379", ]
MA12 <- MA11[!MA11$genes$SystematicName=="TERA03177", ]
MA13 <- MA12[!MA12$genes$SystematicName=="TERA03336", ]
MA14 <- MA13[!MA13$genes$SystematicName=="TERA03412", ]
MA15 <- MA14[!MA14$genes$SystematicName=="TERA03863", ]
MA16 <- MA15[!MA15$genes$SystematicName=="TERA04169", ]

#phase4'''
#w/n array replicate spots
#multiple probe sets mapping to the same gene
#duplicate correlation approach (p94 limma userguide. p10 Smyth 2005 limma)
corfit <- duplicateCorrelation(MA, ndups= 2, spacing="topbottom")
#value should be large & positive
corfit$consensus.correlation

#phase5
#if the M values have different spreads (boxplot) i should normalize between arrays
#limma userguide p76
#Scale the log ratios to have the same median absolute deviation across arrays
#A and M are normalized
MA.Aq <- normalizeBetweenArrays(MA2, method="scale")
plotDensities(MA.Aq)
boxplot(MA.Aq$M~col(MA.Aq$M),names=colnames(MA.Aq$M))
plotMA(MA.Aq)
#Green channel w/ the same empirical distribution across arrays 
#M-ratios are unchanged
MA.Gq <- normalizeBetweenArrays(MA, method="Gquantile")
plotDensities(MA.Gq)
boxplot(MA.Gq$M~col(MA.Gq$M),names=colnames(MA.Gq$M))
plotMA(MA.Gq)

#give the number of spots and number of arrays
dim(MA.Aq)


#phase6	
#bloxplot print out of MA.Aq$M
# P80 < limma Userguide
cols <- MA.Aq$targets$Cy5
cols[cols=="O1"] <- "blue"
cols[cols=="O2"] <- "blue"
cols[cols=="O3"] <- "blue"
cols[cols=="T1"] <- "red"
cols[cols=="T2"] <- "red"
cols[cols=="T3"] <- "red"
cols[cols=="V1"] <- "green"
cols[cols=="V2"] <- "green"
cols[cols=="V3"] <- "green"
cols[cols=="V1c"] <- "green"
cols[cols=="V2c"] <- "green"
cols[cols=="V3c"] <- "green"
cols[cols=="P1"] <- "yellow"
cols[cols=="P2"] <- "yellow"
cols[cols=="P3"] <- "yellow"
cols[cols=="P1c"] <- "yellow"
cols[cols=="P2c"] <- "yellow"
cols[cols=="P3c"] <- "yellow"
cols[cols=="J1"] <- "white"
cols[cols=="J2"] <- "white"
cols[cols=="J1c"] <- "white"
cols[cols=="J2c"] <- "white"
cols[cols=="JA"] <- "orange"
boxplot(MA.Aq$M~col(MA.Aq$M),names=rownames(MA.Aq$targets),col=cols,xlab="Ontogen?se",ylab="M-values")


#phase7
#Simple Comparisons, replicate arrays (Direct two-color experimental design) -- working
fit <- lmFit(MA.Aq)
fit2 <- eBayes(fit)

#phase8
volcanoplot(fit2, highlight=15)

#phase9
#P61 < limma Userguide
topTable(fit2, number=30, adjust="BH")

#===========================

#create a linear matrix using Pool as common RNA source
#setting a design matrix
design <- modelMatrix(targets, ref="Pool")
colnames(design) = c("O1vsPool", "O2vsPool", "O3vsPool",
"T1vsPool", "T2vsPool", "T3vsPool", "V1vsPool", "V2vsPool", "V3vsPool",
"V1cvsPool", "V2cvsPool", "V3cvsPool", "P1vsPool", "P2vsPool", "P3vsPool",
"P1cvsPool", "P2cvsPool", "P3cvsPool", "J1vsPool", "J2vsPool", "J1cvsPool", 
"J2cvsPool", "JAvsPool")
design
fit <- lmFit(MA.Aq, design)


#(trial)linear matrix
contrast.matrix <- makeContrasts(O1, O2, O3, T1, T2, T3, V1, V2, V3, V1c, V2c, V3c, P1, P2, P3, P1c, P2c, P3c, J1, J2, J1c, J2c, JA, 
"O1-T1", "O1-T2", "O1-T3", "O1-V1", "O1-V2", "O1-V3", "O1-V1c", "O1-V2c", "O1-V3c", "O1-P1", "O1-P2", "O1-P3", "O1-P1c", "O1-P2c", "O1-P3c", "O1-J1", "O1-J2", "O1-J1c", "O1-J2c", "O1-JA", 
"O2-T1", "O2-T2", "O2-T3", "O2-V1", "O2-V2", "O2-V3", "O2-V1c", "O2-V2c", "O2-V3c", "O2-P1", "O2-P2", "O2-P3", "O2-P1c", "O2-P2c", "O2-P3c", "O2-J1", "O2-J2", "O2-J1c", "O2-J2c", "O2-JA", 
"O3-T1", "O3-T2", "O3-T3", "O3-V1", "O3-V2", "O3-V3", "O3-V1c", "O3-V2c", "O3-V3c", "O3-P1", "O3-P2", "O3-P3", "O3-P1c", "O3-P2c", "O3-P3c", "O3-J1", "O3-J2", "O3-J1c", "O3-J2c", "O3-JA", 
"T1-V1", "T1-V2", "T1-V3", "T1-V1c", "T1-V2c", "T1-V3c", "T1-P1", "T1-P2", "T1-P3", "T1-P1c", "T1-P2c", "T1-P3c", "T1-J1", "T1-J2", "T1-J1c", "T1-J2c", "T1-JA", 
"T2-V1", "T2-V2", "T2-V3", "T2-V1c", "T2-V2c", "T2-V3c", "T2-P1", "T2-P2", "T2-P3", "T2-P1c", "T2-P2c", "T2-P3c", "T2-J1", "T2-J2", "T2-J1c", "T2-J2c", "T2-JA", 
"T3-V1", "T3-V2", "T3-V3", "T3-V1c", "T3-V2c", "T3-V3c", "T3-P1", "T3-P2", "T3-P3", "T3-P1c", "T3-P2c", "T3-P3c", "T3-J1", "T3-J2", "T3-J1c", "T3-J2c", "T3-JA", 
"V1-V1c", "V1-V2c", "V1-V3c", "V1-P1", "V1-P2", "V1-P3", "V1-P1c", "V1-P2c", "V1-P3c", "V1-J1", "V1-J2", "V1-J1c", "V1-J2c", "V1-JA", 
"V2-V1c", "V2-V2c", "V2-V3c", "V2-P1", "V2-P2", "V2-P3", "V2-P1c", "V2-P2c", "V2-P3c", "V2-J1", "V2-J2", "V2-J1c", "V2-J2c", "V2-JA", 
"V3-V1c", "V3-V2c", "V3-V3c", "V3-P1", "V3-P2", "V3-P3", "V3-P1c", "V3-P2c", "V3-P3c", "V3-J1", "V3-J2", "V3-J1c", "V3-J2c", "V3-JA", 
"V1c-P1", "V1c-P2", "V1c-P3", "V1c-P1c", "V1c-P2c", "V1c-P3c", "V1c-J1", "V1c-J2", "V1c-J1", "V1c-J2c", "V1c-JA", 
"V2c-P1", "V2c-P2", "V2c-P3", "V2c-P1c", "V2c-P2c", "V2c-P3c", "V2c-J1", "V2c-J2", "V2c-J1c", "V2c-J2c", "V2c-JA", 
"V3c-P1", "V3c-P2", "V3c-P3", "V3c-P1c", "V3c-P2c", "V3c-P3c", "V3c-J1", "V3c-J2", "V3c-J1c", "V3c-J2c", "V3c-JA", 
"P1-P1c", "P1-P2c", "P1-P3c", "P1-J1", "P1-J2", "P1-J1c", "P1-J2c", "P1-JA", 
"P2-P1c", "P2-P2c", "P2-P3c", "P2-J1", "P2-J2", "P2-J1c", "P2-J2c", "P2-JA", 
"P3-P1c", "P3-P2c", "P3-P3c", "P3-J1", "P3-J2", "P3-J1c", "P3-J2c", "P3-JA", 
"P1c-J1", "P1c-J2", "P1c-J1c", "P1c-J2c", "P1c-JA", 
"P2c-J1", "P2c-J2", "P2c-J1c", "P2c-J2c", "P2c-JA", 
"P3c-J1", "P3c-J2", "P3c-J1c", "P3c-J2c", "P3c-JA", 
"J1-J1c", "J1-J2c", "J1-JA", 
"J2-J1c", "J2-J2c", "J2-JA", 
"J1c-JA", 
"J2c-JA", levels=design)
contrast.matrix


#(trial)works with the matrix above (gives a degrees of freedom error) -- not resolved yet
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#the ordinary t-statistics can be recoverd after fitting a linear model
tstat.ord <- fit$coef/fit$stdev.unscaled/fit$sigma

#spatial image plot
#bioC userguide
targets <- rnorm(9*3*18*18)
imageplot(targets,layout=list(ngrid.r=9,ngrid.c=3,nspot.r=18,nspot.c=18))

	



