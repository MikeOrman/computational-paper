#--------------------------FUNCTION FOR WRITING CBIO SAMPLES--------------------
#dataset: Col1 = Hugo_Symbols
#info[1] = study, info[2] = output filename
cBio_SampleFile <- function(dataset, info){
  names = colnames(dataset)[2:ncol(dataset)]
  format = c()
  format[1:length(names)] = ":"
  study = c()
  study[1:(ncol(dataset)-1)] = info[1]
  sample = data.frame(study, format, names)
  write.table(sample, file = info[2], 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep="")
  return(sample)
}