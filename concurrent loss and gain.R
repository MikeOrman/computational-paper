source("generank.R")
concurrent.loss.data <- function(subset, dataset){
  #Calc exclusive copy number loss genes
  subset.loss = lossRank(subset)[,c(1,3,4,5)]
  colnames(subset.loss) <- c("Hugo_Symbol", "Subset loss", "Subset unaltered", "Subset loss ratio")
  dataset.loss = lossRank(dataset)[,c(1,3,4,5)]
  colnames(dataset.loss) <- c("Hugo_Symbol", "Dataset loss", "Dataset unaltered", "Dataset loss ratio")
  both.loss = merge.data.frame(subset.loss, dataset.loss)
  for (i in 1:nrow(both.loss)){
    normalized.loss = both.loss$`Subset loss ratio`[i] / both.loss$`Dataset loss ratio`[i]
    matrix = matrix(as.numeric(both.loss[i, c(2, 5, 3, 6)]), ncol = 2)
    test = fisher.test(matrix)$p.value
    both.loss[i,8] = normalized.loss
    both.loss[i,9] = test
  }
  both.loss[,10] <- p.adjust(both.loss[,9], method = "fdr")
  colnames(both.loss)[10] <- "adjusted loss p val"
  colnames(both.loss)[8] <- "normalized subtype loss ratio"
  colnames(both.loss)[9] <- "loss p val"
  return(both.loss)
}

concurrent.gain.data <- function(subset, dataset){
  #Calc exclusive copy number gain genes  
  subset.gain = gainRank(subset)[,c(1,3,4,5)]
  colnames(subset.gain) <- c("Hugo_Symbol", "Subset gain", "Subset unaltered", "Subset gain ratio")
  dataset.gain = gainRank(dataset)[,c(1,3,4,5)]
  colnames(dataset.gain) <- c("Hugo_Symbol", "Dataset gain", "Dataset unaltered", "Dataset gain ratio")
  both.gain = merge.data.frame(subset.gain, dataset.gain)
  for (i in 1:nrow(both.gain)){
    normalized.gain = both.gain$`Subset gain ratio`[i] / both.gain$`Dataset gain ratio`[i]
    matrix = matrix(as.numeric(both.gain[i, c(2, 5, 3, 6)]), ncol = 2)
    test = fisher.test(matrix)$p.value
    both.gain[i,8] = normalized.gain
    both.gain[i,9] = test
  }
  both.gain[,9] <- p.adjust(both.gain[,9], method = "fdr")
  colnames(both.gain)[9] <- "adjusted gain p val"
  colnames(both.gain)[8] <- "normalized subtype gain ratio"
  return(both.gain)
}