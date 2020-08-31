library(tidyverse)

# patients with pre and post cancer cells
patient_cells = read.table("FELINE_cancer_cells.txt", sep="\t", header=T)
patient_cells %>% filter(Pre %in% c("Yes"), Post %in% c("Yes")) -> patient_cells
patient_cells %>% filter(!PID %in% c("P04", "P32")) -> patient_cells

# clinical sorting
clinic = read.table("FELINE_patient_1_46.clinical_data.txt", sep="\t", header=T)
clinic %>% select(Sample, orig.ident, ARM, Response) -> clinic
clinic$PID = gsub("FEL0","P", clinic$orig.ident)

# merge and sort 
clinic_sort = merge(clinic, patient_cells, on=PID)
clinic_sort %>% arrange(Response, ARM, PID) -> clinic_sort
print(clinic_sort)
write.table(clinic_sort, "FELINE_patient_info.sorted.txt", quote=F, row.names=F, col.names=T, sep="\t")
