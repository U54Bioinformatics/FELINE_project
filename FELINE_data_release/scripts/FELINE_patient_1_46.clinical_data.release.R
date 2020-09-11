library(tidyverse)

x = read.table("FELINE_patient_1_46.clinical_data.txt", header=T, sep="\t")
x$Patient.ID = gsub('FEL0', 'P', x$orig.ident)
x %>% select(Patient.ID, orig.ident, ARM, Treatment, Response) %>% filter(!Patient.ID %in% c("P18")) -> x
write.table(x, "FELINE_patient_1_46.clinical_data.release.txt", row.names=F, col.names=T, quote=F, sep="\t")

