
## Launch 01funcs.R

#########################################################################
#	SECTION -1-								#
#########################################################################
ratio <- matrix(sample(seq(-5, 5, by=0.1), 100, replace=T), 100, 4, dimnames=list(paste("g", 1:100, sep=""), paste("DEG", 1:4, sep="")), byrow=T)
   # Creates a sample matrix of gene expression log2 ratios. This could be any data type!
setlistup <- sapply(colnames(ratio), function(x) rownames(ratio[ratio[,x]>=1,]))
setlistdown <- sapply(colnames(ratio), function(x) rownames(ratio[ratio[,x]<=-1,]))
   # Identifies all genes with at least a two fold up or down regulation and stores the corresponding gene identifiers
   # in setlistup and setlistdown, respectively.
OLlistup <- overLapper(setlist=setlistup, sep="_", type="vennsets")
OLlistdown <- overLapper(setlist=setlistdown, sep="_", type="vennsets")
counts <- list(sapply(OLlistup$Venn_List, length), sapply(OLlistdown$Venn_List, length))
vennPlot(counts=counts, ccol=c("red", "blue"), colmode=2, mysub="Top: DEG UP; Bottom: DEG Down", yoffset=c(0.3, -0.2))


source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R") # Imports required functions.
setlist <- list(A=sample(letters, 18), B=sample(letters, 16), C=sample(letters, 20), D=sample(letters, 22), E=sample(letters, 18), F=sample(letters, 22, replace=T))
setlist5 <- setlist[1:5]
OLlist5 <- overLapper(setlist=setlist5, sep="_", type="vennsets")
counts <- sapply(OLlist5$Venn_List, length)
vennPlot(counts=counts, ccol=c(rep(1,30),2), lcex=1.5, ccex=c(rep(1.5,5), rep(0.6,25),1.5)) # Plots a non-proportional 5-way Venn diagram.
