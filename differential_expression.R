library(limma)
differential_expression <- function(eset){
mRNA_data = as.data.frame(exprs(eset))
sample_data = pData(eset)
feature_data = fData(eset)
design = model.matrix(~0 + as.factor(subtype),
                      data = sample_data)
if ((sum(eset@phenoData@data[["subtype"]] == "subtype1") > 0) &
    (sum(eset@phenoData@data[["subtype"]] == "subtype2") > 0) &
    (sum(eset@phenoData@data[["subtype"]] == "unclassified") > 0)){
  colnames(design) <- c("subtype1", "subtype2", "unclassified")
  } else {colnames(design) <- c("subtype1", "subtype2")}
contrast = makeContrasts(subtype1vssubtype2 = subtype1-subtype2, levels = design)
#Calculated moderated statistics
fit = lmFit(mRNA_data, design)
contrasts = contrasts.fit(fit, contrast)
fit.statistics = eBayes(contrasts, trend = TRUE)
sigtable = topTable(fit.statistics, adjust = 'fdr', number = nrow(mRNA_data),
                    confint = T)
sigtable[,(ncol(sigtable)+1)] = rownames(sigtable)
subtype1_count = sum(eset@phenoData@data$subtype == "subtype1")
subtype2_count = sum(eset@phenoData@data$subtype == "subtype2")
output = data.frame(sigtable$V9, subtype1_count, subtype2_count, 
                    sigtable$adj.P.Val, sigtable$logFC)
colnames(output) = c("Hugo_Symbol", "Subtype 1 mRNA sample count", "Subtype 2 mRNA sample count",
                     "Moderated FDR", "log2FC")
return(output)
}