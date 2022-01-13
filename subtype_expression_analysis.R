#-------------DETERMINE EXPRESSION FC/SIGNIFCANCE OF LOST ENRICHED GENES------------
source("mRNA_eset.R")
source("differential_expression.R")
source("subtype_subset.R")
subtype_with_loss_expression <- function (subtype, loss_genes, mRNA, comparison){
  output = data.frame()
  for (i in 1:length(loss_genes)){
    print(i)
    subtype_with_loss = subtype_subset(subtype, c(loss_genes[i], comparison[1]))
    subtype_without_loss = subtype_subset(subtype, c(loss_genes[i], comparison[2]))
    if (class(subtype_with_loss) == "data.frame" & class(subtype_without_loss) == "data.frame"){
      genes = subtype$Hugo_Symbol
      names1 = colnames(subtype_with_loss)[2:ncol(subtype_with_loss)]
      names2 = colnames(subtype_without_loss)[2:ncol(subtype_without_loss)]
      eset = mRNA_eset(mRNA, subtype_with_loss, subtype_without_loss)
      if (class(eset) == "data.frame"){
        output[i,1] = loss_genes[i]
        output[i,2] = eset[1,1]
        output[i,3] = eset[1,2]
        output[i,4] = NA
        output[i,5] = NA
        } else {
          exp = differential_expression(eset)
          index = which(loss_genes[i] == exp$Hugo_Symbol)
          gene_exp = exp[index,]
          output[i,1] = gene_exp[1,1]
          output[i,2] = gene_exp[1,2]
          output[i,3] = gene_exp[1,3]
          output[i,4] = gene_exp[1,4]
          output[i,5] = gene_exp[1,5]
        }
    }
    if (class(subtype_with_loss) == "character" & class(subtype_without_loss) == "data.frame"){
      output[i,1] = loss_genes[i]
      output[i,2] = 0
      output[i,3] = ncol(subtype_without_loss) - 1
      output[i,4] = NA
      output[i,5] = NA
      }
    if (class(subtype_with_loss) == "data.frame" & class(subtype_without_loss) == "character"){
      output[i,1] = loss_genes[i]
      output[i,2] = ncol(subtype_with_loss) - 1
      output[i,3] = 0
      output[i,4] = NA
      output[i,5] = NA
    }
  }
    colnames(output) <- c("Hugo_Symbol", paste(comparison[1], "mRNA sample count"),
                          paste(comparison[2],"mRNA sample count"),
                          "Adjusted p-value", "logFC")
    return(output)
}