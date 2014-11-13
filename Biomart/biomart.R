## from (durinck2009mapping)
library("biomaRt")
listMarts()
{
## choose a DATABASE
ensembl <- useMart("ensembl")
## choose a DATASET
listDatasets(ensembl)
## Select a DATASET
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
}
## FILTERS of the selected DATASET
filters <- listFilters(ensembl)
filters
## ATTRIBUTES
attributes <- listAttributes(ensembl)
attributes

customIDs=c("OPHN3") ## protein ID (from uniprot)
{
getBM(attributes=c('hgnc_symbol', 'go_id','external_gene_id','description'), 					
				filters = c( "hgnc_symbol","with_protein_id"), 
				values = list(customIDs,TRUE), 
				mart=ensembl)
}

 



























