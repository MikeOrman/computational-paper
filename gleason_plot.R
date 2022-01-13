#dataset 1 = loss subtype 
#dataset 2 = wt/gain subtype
#gleason data = dataframe of sample_ID, gleason score
#dataset.names = character vector of length 2 for dataset 1 and dataset 2 genotypes
source("gleason.R")
gleason_bar_chart <- function(dataset1, dataset2, gleason.data, dataset.names){
  dataset1.gleason = gleason_data(dataset1, gleason.data)$`Stratified by risk group`
  dataset1.gleason[,4] = dataset.names[1]
  dataset2.gleason = gleason_data(dataset2, gleason.data)$`Stratified by risk group`
  dataset2.gleason[,4] = dataset.names[2]
  combined = rbind(dataset1.gleason, dataset2.gleason)
  combined = na.omit(combined)
  colnames(combined)[4] <- "dataset"
  return(combined)
}