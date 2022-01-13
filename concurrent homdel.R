source("homdelRank.R")
concurrent.homdel.data <- function(subset, dataset){
  #Calc exclusive copy number loss genes
  subset.loss = homdelRank(subset)[,c(1,3,4,5)]
  colnames(subset.loss) <- c("Hugo_Symbol", "Subset homdel", "Subset unaltered", "Subset homdel ratio")
  dataset.loss = lossRank(dataset)[,c(1,3,4,5)]
  colnames(dataset.loss) <- c("Hugo_Symbol", "Dataset homdel", "Dataset unaltered", "Dataset homdel ratio")
  both.loss = merge.data.frame(subset.loss, dataset.loss)
  for (i in 1:nrow(both.loss)){
    matrix = matrix(as.numeric(both.loss[i, c(2, 5, 3, 6)]), ncol = 2)
    test = fisher.test(matrix)$p.value
    both.loss[i,8] = test
  }
  both.loss[,8] <- p.adjust(both.loss[,8], method = "fdr")
  colnames(both.loss)[8] <- "adjusted loss p val"
  return(both.loss)
}