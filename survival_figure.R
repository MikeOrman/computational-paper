d = read.table("data_clinical_patient.txt", sep = "\t")
Patient_ID = as.numeric(d[1:128, 1])
Month = as.numeric(d[1:128, 6])
Death = d[1:128, 5]
Death[Death == "0:LIVING"] = 0
Death[Death == "1:DECEASED"] = 1
Death = as.numeric(Death)
patient_data = data.frame(Patient_ID, Death, Month)
sample_data = read.table("abida.sample.txt", header = TRUE, sep = "\t")
IDs = sample_data[,c(1,2)]
source("survival.R")
source("subtype_subset")
dataset1 = subtype_subset(met.coloss, c("SUFU", "loss"))
dataset2 = subtype_subset(met.coloss, c("SUFU", "WT/gain"))
subtype = c("SUFU loss", "SUFU WT/gain")
surv_data = surv_analysis(dataset1, dataset2, patient_data, IDs)
plot(surv_data$fit, ylab = "Overall Survival", xlab = "Months", main = "Survival Analysis",
     col = c("Blue", "Dark Red"), lwd = 2.5, 
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
legend("topright", 
       legend = c(subtype[1], subtype[2], sprintf("p = %s", round(surv_data$sd$chisq, 4))),
       col = c("Blue", "Dark Red", "Black"), lty=1:1, cex=1, bty = "n")
