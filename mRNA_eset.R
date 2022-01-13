#Creates an expression set of mRNA transcriptomic data for two groups from an experiment
#mRNA = mRNA data of the experiment
#first subtype dataset
#second subtype dataset
library(Biobase)
mRNA_eset <- function (mRNA, dataset1, dataset2){
  names1 = colnames(dataset1)[2:ncol(dataset1)]
  names2 = colnames(dataset2)[2:ncol(dataset2)]
  genes = dataset1$Hugo_Symbol
  mRNA = mRNA[(duplicated(mRNA$Hugo_Symbol) == FALSE),]
  mRNA = na.omit(mRNA)
  index = c()
  for (i in 1:nrow(mRNA)){
    if (sum(mRNA$Hugo_Symbol[i] == genes) == 1) {
      index[i] = TRUE
    } else {index[i] = FALSE}
  }
  mRNA <- mRNA[index,]
  #Subsets mRNA data by provided group names
  index1 = c(TRUE)
  for (i in 2:ncol(mRNA)) {
    if (sum(colnames(mRNA)[i] == names1) > 0) {
      index1[i] = TRUE
    }
    if ((sum(colnames(mRNA)[i] == names1) > 0) == FALSE) {
      index1[i] = FALSE
    }
  }
  mRNA1 = mRNA[, index1]
  index2 = c(TRUE)
  for (i in 2:ncol(mRNA)) {
    if (sum(colnames(mRNA)[i] == names2) > 0) {
      index2[i] = TRUE
    }
    if ((sum(colnames(mRNA)[i] == names2) > 0) == FALSE) {
      index2[i] = FALSE
    }
  }
  mRNA2 = mRNA[, index2]
  #mRNA samples not in either group
  if (class(mRNA1) == "data.frame" & class(mRNA2) == "data.frame") {
    group.names=c(colnames(mRNA1)[2:ncol(mRNA1)], colnames(mRNA2)[2:ncol(mRNA2)])
  }
  if (class(mRNA1) == "data.frame" & class(mRNA2) == "character") {
    group.names=colnames(mRNA1)[2:ncol(mRNA1)]
  }
  if (class(mRNA1) == "character" & class(mRNA2) == "data.frame") {
    group.names=colnames(mRNA2)[2:ncol(mRNA2)]
  }
  if (class(mRNA1) == "character" & class(mRNA2) == "character") {
    group.names=NA
  }
  index3 = c(TRUE)
  for (i in 2:ncol(mRNA)){
    if (sum(colnames(mRNA)[i] == group.names) == 0) {
      index3[i] = TRUE
    } else {index3[i] = FALSE}
  }
  mRNA3 = mRNA[, index3]
  #Check if mRNA subtype datasets contain at least 3 samples
  if (sum(index1) >= 4 & sum(index2) >=4){
    #Read Sample Data
    if (sum(ncol(mRNA1) > 1) > 0) {
      subtype1 = c()
      subtype1[1:(ncol(mRNA1)-1)] = "subtype1"
    }
    if (sum(ncol(mRNA2) > 1) > 0){
      subtype2 = c()
      subtype2[1:(ncol(mRNA2)-1)] = "subtype2"
    }
    if (sum(ncol(mRNA3) > 1) > 0){
      subtype3 = c()
      subtype3[1:(ncol(mRNA3)-1)] = "unclassified"
    }
    if (exists("subtype1") & exists("subtype2") & exists("subtype3")) {
      exprs = merge.data.frame(mRNA1, mRNA2)
      exprs = merge.data.frame(exprs, mRNA3)
      rownames(exprs) = mRNA$Hugo_Symbol
      exprs = exprs[,2:ncol(exprs)]
      exprs = as.matrix(exprs)
      subtype = c(subtype1, subtype2, subtype3)
      pData = data.frame(subtype)
      rownames(pData)[1:(ncol(mRNA1)-1)] <- colnames(mRNA1)[2:ncol(mRNA1)]
      rownames(pData)[ncol(mRNA1):((ncol(mRNA1)+ncol(mRNA2))-2)] <- colnames(mRNA2)[2:ncol(mRNA2)]
      rownames(pData)[((ncol(mRNA1)+ncol(mRNA2))-1):(((ncol(mRNA1)+ncol(mRNA2))-3)+ncol(mRNA3))] <- colnames(mRNA3)[2:ncol(mRNA3)]
      pData = AnnotatedDataFrame(pData)
    }
    if (exists("subtype1") & exists("subtype2") & (exists("subtype3") == FALSE)) {
      exprs = merge.data.frame(mRNA1, mRNA2)
      rownames(exprs) = mRNA$Hugo_Symbol
      exprs = exprs[,2:ncol(exprs)]
      exprs = as.matrix(exprs)
      subtype = c(subtype1, subtype2)
      pData = data.frame(subtype)
      rownames(pData)[1:(ncol(mRNA1)-1)] <- colnames(mRNA1)[2:ncol(mRNA1)]
      rownames(pData)[ncol(mRNA1):((ncol(mRNA1)+ncol(mRNA2))-2)] <- colnames(mRNA2)[2:ncol(mRNA2)]
      pData = AnnotatedDataFrame(pData)
    }
    #Formats feature data in raw data variable
    fData = data.frame(rownames(exprs),
                       row.names = rownames(exprs))
    colnames(fData) = "Hugo_Symbol"
    fData = AnnotatedDataFrame(fData)
    #Creates and saves expression set of group 1 and group 2 mRNA data
    mRNA.eset <- ExpressionSet(assayData=exprs, phenoData = pData, featureData = fData)
    return(mRNA.eset)
  }
  #Returns sample sizes if either mRNA subtype group has less than 3 samples
  else {
    output = data.frame()
    output[1,1] = (sum(index1)-1)
    output[1,2] = (sum(index2)-1)
    colnames(output) = c("Samples subtype 1", "Samples subtype 2")
    return(output)
  }
}