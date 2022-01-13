#sample.names = vector of sample names
#gleason.df is a  dataframe where, 
#         column 1 = sample names of all samples in a study
#         column 2 = corresponding gleason scores for all samples the study
library("dplyr")
gleason_data <- function(samples, gleason.df){
  sample.names = colnames(samples)[2:ncol(samples)]
  index = c()
  for (i in 1:nrow(gleason.df)) {
    if ((sum(which(sample.names[i] == gleason.df))) > 0) {
    index[i] = which(sample.names[i] == gleason.df)
    }
  }
  index = na.omit(index)
  output = data.frame(gleason.df[index, 1], gleason.df[index, 2])
  colnames(output) <- c("Sample_ID", "Gleason Score")
  output <- output  %>%
    group_by(`Gleason Score`) %>%
    dplyr::summarize(n())
  output = as.data.frame(output)
  output[,1] = as.numeric(output[,1])
  output = output
  total = sum(output$`n()`)
  percentage = output[,2]/total
  output[,3] = percentage
  colnames(output)[3] <- "Percentage"
  output = output[order(output$`Gleason Score`, decreasing = FALSE),]
  risk = c("Low", "Medium", "High")
  count = c()
  if (sum((output$`Gleason Score`== 6) == TRUE) > 0) {
    count[1] = output[output$`Gleason Score`== 6,2]
  } else {count[1] = 0}
  if (sum((output$`Gleason Score`== 7) == TRUE) > 0) {
    count[2] = output[output$`Gleason Score`== 7,2]
  } else {count[2] = 0}
  if (sum((output$`Gleason Score`== 8) == TRUE) > 0) {
    count[3] = output[output$`Gleason Score`== 8,2]
  } else {count[3] = 0}
  if (sum((output$`Gleason Score`== 9) == TRUE) > 0) {
    count[3] = count[3] + output[output$`Gleason Score`== 9,2]
  }
  percentage = count / sum(count)
  output2 = data.frame(risk, count, percentage)
  output3 = list(output, output2)
  names(output3) <- c("Stratified by score", "Stratified by risk group")
  return(output3)
}