#Col 1 must be named "Hugo_Symbol"
library("AnnotationDbi")
library("org.Hs.eg.db")
cancer_annotations <- function(dataset){
  cols = ncol(dataset)
  oncogenes = read.table("cancer_gene_census.csv", sep = ",", header = TRUE)
  for (i in 1:nrow(dataset)){
    if (sum(dataset[i,1] == oncogenes[,1], na.rm = TRUE) > 0){
      dataset[i,(cols+1)] = TRUE
    }
    else {
      dataset[i,(cols+1)] = FALSE
    }
  }
  TSGs = read.table("Human_TSGs.txt", sep = "\t", header = TRUE)
  for (i in 1:nrow(dataset)){
    if (sum(dataset[i,1] == TSGs[,2], na.rm = TRUE) > 0){
      dataset[i,(cols+2)] = TRUE
    }
    else {
      dataset[i,(cols+2)] = FALSE
    }
  }
  colnames(dataset)[(ncol(dataset)-1)] <- "Oncogene"
  colnames(dataset)[(ncol(dataset))] <- "TSG"
  return(dataset)
}
loci <- function(dataset){
  uniKeys <- keys(org.Hs.eg.db, keytype="UNIPROT")
  cols <- c("SYMBOL", "MAP")
  band = select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="UNIPROT")[,2:3]
  columns = ncol(dataset)
  colnames(band) = c("Hugo_Symbol", "Band")
  for (i in 1:nrow(dataset)){
    index = which(dataset$Hugo_Symbol[i] == band$Hugo_Symbol)
    index = index[1]
    dataset[i, columns + 1] = band$Band[index]
  }
  colnames(dataset)[ncol(dataset)] <- "Band"
  return(dataset)
}
