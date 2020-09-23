library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(cowplot)

# patient id
patient_id = read.table("FELINE_patient_1_46.clinical_data.with_three_IDs.txt", sep="\t", header=T)
patient_id %>% select(Patient.ID, orig.ident, Patient.Study.ID) -> patient_id
names(patient_id) <- c("patient", "orig.ident", "Patient.Study.ID")

t = read.table("FELINE_gene_cnv_barplot.patient.STable_allCNV.txt", sep="\t", header=T)
t_m = merge(t, patient_id, on=patient)
t_m %>% arrange(variable, patient) -> t_m
write.table(t_m, "FELINE_gene_cnv_barplot.patient.STable_allCNV.three_IDs.txt", sep="\t", quote=F, row.names=F)
