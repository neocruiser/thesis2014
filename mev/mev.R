

## Processing output data of MEV for microarrays (the follwing code belong to preprocessing file)

# ========================================================================================
# Briefing before processing MeV output files
# ========================================================================================
# (step 1) = paste only CUSTs (MeV) into notepad++. Save File1
# (step 2) = type I file = retrieve Locus ID of CUST. Save File2
# (step 3) = type II file = retrieve nucleotide sequence for significant MeV spots. Save File3
# (step 4) = type III file = remove sequence, keep locus|GO ID|name. Save File4
# (additional step 5) = paste ANOVA data (MeV). Save File0
# (additional step 6) = import File0 and File2 in excel


## START
# ========================================================================================
# Reason: create a type I file. 
# ========================================================================================
# Retrieve "locus names" to "CUSTs" selected with MeV
# (1) add "cust" column in cust.txt files

library(limma)
setwd("C:/Dropbox/Workshop2013/Work/Microarrays/10- MeV/Gene Selection/")
id.t <- readTargets("id.txt")         #from bacterial spot removal
# (2) change only file names then repeat the following codes
CUST <- readTargets("anova1.txt")   #input file (name); CUSTs generated w/ MeV
CUST.t <- CUST[!duplicated(CUST$cust),]    #remove duplicates, hence creating a data frame of "characters" which could be read by %in%
locus <- id.t[id.t$cust %in% c(CUST.t),]
x <- as.data.frame(locus)             #convert to "data frame" because "%in% locus" is a "character matrix"
locus.t <- x[!duplicated(x$locus),]         #only data frames could be treated for duplication
write.table(locus.t, "anova2.txt", sep="\t")
rm(CUST, CUST.t, locus, x, locus.t)
## END


## START
# ========================================================================================
# Reason: Create type II files. w/o SAM. w/o CUST. Retrieve "Nucleotide sequence"
# ========================================================================================
# Find nucleotide sequences from differentially expressed genes (CUSTs from MeV output)
# edit the "export_seq_blast15457.fasta" (build a data frame that can be imported in R, reason to find sequences specific to MeV CUSTs)
# (0) check that there is only 15457 occurences being replaced
# (1) notepad++ search: (^>.*$)\r\n([^>]*$) and replace: \1\t\2 (in regular expression mode)
# (2) notepad++ search: \| and replace: \t\|
# (3) notepad++ search: ^> and replace: >\t
# (4) add first line (x locus description sequence)
# (5) save as seperate file "exportID.txt"
# (6) when finished, reverse repeat steps (3), (2), (1)
# (7) notepad++ search: " and replace: 'blank'
# (8) notepad++ search: ^>\t and replace: ^>
# (9) notepad++ search: \t\| and replace: \|
# (10) notepad++ search: (^>.*)\t and replace: \1\r\n
# (hint) Procedure done once per dataset. Meaning exportID is for all sequences

library(limma)
setwd("C:/Dropbox/Workshop2013/Work/Microarrays/10- MeV/Gene Selection/")
ID <- readTargets("exportID.txt")        # containing the 15k annotated Blast2GO sequences
# (7) change only file names then repeat the following codes
exp <- read.table("anova2.txt", header=TRUE, sep="")      # containing the MeV exported experimental locuses (SAM statistics)
locus <- as.character(exp[,2])
seq <- ID[ID$locus %in% c(locus),]
write.table(seq, "anova2.txt", sep="\t")
rm(exp, seq, locus)
## END

## START
# ========================================================================================
# Reason: create type III files
# ========================================================================================
# Traitement des séquences après les avoirs indexées puis extraits, create a tab delimited files
# (O) copy AnnotSeq.txt files into another directory
# (update 1) notepad++ search: .*?> and replace: >
# (update 2) notepad++ search: " and replace: blank
# (update 3) notepad++ search: >\t and replace: >
# (update 4) notepad++ search: \t| and replace: |
# (update 5) notepad++ search: \t and replace: \r\n
# (1) notepad++ search: ^[^>]*$ and replace: 'blank' (removes all nucleotide sequences)
# (2) notepad++ search: \r\n> and replace: 'blank' (removes all blank lines)
# (4) notepad++ search: | and replace: \t (adds a tab between variables)
# (not necessary anymore) notepad++ search: (?<=\t)(.*$) and replace: "\1"
# (6) notepad++ add header names (oeuf  trocophor)(oeuf veliger)(oeuf pediveliger) etc.
# (hint) All files must be sep=""
# (tip) R will take every word seperated by a space as a new variable then prompt an error, this is why you should add quotation mark to every observation 
# (7) make sure you have:
# (A) File set for CUST only
# (B) File set for CUST vs locus
# (C) File set for locus vs descriptin (w/o sequence)
# (D) File set for SAM MeV data
# (E) File set for CUST - locus - description - sequence - SAM
# (hint) Merge B - C - D into one excel file
# (tip) remove space between 2 words, if so, add quotation marks (eg, "hello world")
# (tip) type I files contains 3 variables, 1 of which should be deleted.
# (1) notepad++ search: ^.*\t(CUST.*?$) and replace: \1
# type I files == CompleteCust.txt
# type II files == CompleteLocus.txt
# type III files == CompleteSAM.txt

library(limma)
setwd("C:/Dropbox/Workshop2013/Work/Microarrays/10- MeV/Gene Selection/")
data.path <- "C:/Dropbox/Workshop2013/Work/Microarrays/10- MeV/Gene Selection/"
datafiles <- list.files(path=".", pattern="locus_.*?.txt$")
data <- do.call(rbind, lapply(datafiles, read.table))
write.table(data, "CompleteCust.txt", sep="\t")
rm(datafiles, data)