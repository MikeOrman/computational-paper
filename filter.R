#Filters genes that are lost at >=10%
#Outputs a data frame that can be merged with full dataset to reduce dimentions
source("generank.R")
loss_filter <- function(subset){
  subset.loss = lossRank(subset)[,c(1,3,4,5)]
  filtered_genes = as.data.frame(subset.loss$Hugo_Symbol[subset.loss$`copy number loss freq` >= .15])
  colnames(filtered_genes) <- "Hugo_Symbol"
  return(filtered_genes)
}