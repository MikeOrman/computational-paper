library(survival)
#dataset1 and dataset 2 header: Col 1 = Hugo Symbol, Col 2 to end = Sample Names
#patient_data header: Col 1 = Patient ID, Col 2 - Event Status, Col 3 = Month
#IDs header: Col 1 = Sample ID, Col 2 = Patient ID
surv_analysis <- function(dataset1, dataset2, patient_data, IDs){
  set1_sample_names = colnames(dataset1)[2:ncol(dataset1)]
  set1_ID = c()
  set1_status = c()
  set1_time = c()
  set1_group = c()
  for (i in 1:length(set1_sample_names)){
    set1_ID[i] = IDs[which(set1_sample_names[i] == IDs[,1]),2]
    if (identical(which(set1_ID[i] == patient_data[,1]), integer(0))){
      set1_status[i] = NA
      set1_time[i] = NA
      set1_group[i] = "Set 1"}
      else {
        set1_status[i] = patient_data[which(set1_ID[i] == patient_data[,1]),2]
        set1_time[i] = patient_data[which(set1_ID[i] == patient_data[,1]),3]
        set1_group[i] = "Set 1"}
  }
  set1_survival_data = data.frame(set1_ID, set1_group, set1_time, set1_status)
  colnames(set1_survival_data) <- c("Patient_ID", "Group", "Time", "Status")
  set2_sample_names = colnames(dataset2)[2:ncol(dataset2)]
  set2_ID = c()
  set2_status = c()
  set2_time = c()
  set2_group = c()
  for (i in 1:length(set2_sample_names)){
    set2_ID[i] = IDs[which(set2_sample_names[i] == IDs[,1]),2]
    if (identical(which(set2_ID[i] == patient_data[,1]), integer(0))){
      set2_status[i] = NA
      set2_time[i] = NA
      set2_group[i] = "Set 2"}
    else {
      set2_status[i] = patient_data[which(set2_ID[i] == patient_data[,1]),2]
      set2_time[i] = patient_data[which(set2_ID[i] == patient_data[,1]),3]
      set2_group[i] = "Set 2"}
  }
  set2_survival_data = data.frame(set2_ID, set2_group, set2_time, set2_status)
  colnames(set2_survival_data) <- c("Patient_ID", "Group", "Time", "Status")
  survival_data = rbind(set1_survival_data, set2_survival_data)
  survival_data = na.omit(survival_data)
  surv_object <- Surv(time = survival_data$Time, event = survival_data$Status)
  fit <- survfit(surv_object ~ Group, data = survival_data)
  sd <- survdiff(surv_object ~ Group, data = survival_data)
  output = list(fit, sd)
  names(output) <- c("fit", "sd")
  return(output)
}