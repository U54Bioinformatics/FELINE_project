library(tidyverse)

# clinic info
clinic = read.table("FELINE_patient_1_46.clinical_data.txt", sep="\t", header=T)
clinic %>% select(Sample, orig.ident, Patient.Study.ID, ARM, Response) -> clinic

# tumor size prediction from Jason
tumor_size = read.table("FELINE_clinical_response.txt", sep="\t", header=T)
tumor_size %>% select(Patient.Study.ID, Day, Burden, Burden_t0) -> tumor_size

# merge
clinic_tumor_size = merge(clinic, tumor_size, on=Patient.Study.ID)
clinic_tumor_size$Burden_change = clinic_tumor_size$Burden/clinic_tumor_size$Burden_t0
write.table(clinic_tumor_size, "FELINE_clinical_response.size_change.txt", quote=F, row.names=F, col.names=T, sep="\t")


