
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(cowplot)

# patient id
patient_id = read.table("FELINE_patient_1_46.clinical_data.with_three_IDs.txt", sep="\t", header=T)
patient_id %>% select(Patient.ID, orig.ident, Patient.Study.ID) -> patient_id
names(patient_id) <- c("patient", "orig.ident", "Patient.Study.ID")

# gene list
gene_list_gain <- c("AKT1", "AKT3", "CCND1", "CCNE2", "CDK6", "ERBB4", "FGFR1", "FGFR2", "MYC", "PIK3CA")
gene_list_gain_label <- paste0(gene_list_gain, "_gain")
gene_list_loss <- c("ESR1", "RB1", "TP53")
gene_list_loss_label <- paste0(gene_list_loss, "_loss")
gene_list <- c(gene_list_gain, gene_list_loss)
gene_list_label <- c(gene_list_gain_label, gene_list_loss_label)

# cnv call
x = read.table("FELINE_FACETS.gene_cnv_call.short_list.STable_CNVcall.txt", header=T, sep="\t")
x$patient <- substr(x$sample, 1, 3)
x %>% select(c("patient", "response", gene_list)) -> x_NR
x_NR
melt(x_NR, id.vars=c("patient", "response")) -> x_NR_m
x_NR_m
## set gain and loss status
x_NR_m$value[x_NR_m$value=="Gain"] = "gain"
x_NR_m$value[x_NR_m$value=="LOH"] = "loss"
x_NR_m$value[x_NR_m$value=="Loss"] = "loss"
## create label
x_NR_m$cnv = paste(x_NR_m$variable, x_NR_m$value, sep="_")
write.table(x_NR_m, "FELINE_gene_cnv_barplot.patient.STable_CNVcall.txt", quote=F, sep="\t")
x_NR_m_distinct = distinct(x_NR_m)
x_NR_m_distinct %>% arrange(x_NR_m_)
write.table(x_NR_m_distinct, "FELINE_gene_cnv_barplot.patient.STable_allCNV.txt", quote=F, sep="\t", row.names=F)
#x_NR_m
x_NR_m %>% filter(cnv %in% gene_list_label) -> x_NR_m_data
x_NR_m_data <- distinct(x_NR_m_data)
x_NR_m_data
write.table(x_NR_m_data, "FELINE_gene_cnv_barplot.patient.STable_targetCNV.txt", quote=F, sep="\t", row.names=F)
x_NR_m_data_table <- data.frame(table(x_NR_m_data$response, x_NR_m_data$cnv))
names(x_NR_m_data_table) <- c("Response", "CNV", "Target_mutation")
x_NR_m_data_table$Others <- 12 - x_NR_m_data_table$Target_mutation
x_NR_m_data_table
#x_NR_m_data_table$Label <- paste0(x_NR_m_data_table$CNV, " (", as.integer(100*(x_NR_m_data_table$Num/sum(x_NR_m_data_table$Num))), "%)")
#x_NR_m_data_table

# prepare fisher test df
x_NR_m_data_table %>% filter(Response %in% c("Responder")) -> x_NR_m_data_table_R
x_NR_m_data_table %>% filter(Response %in% c("Non-responder")) -> x_NR_m_data_table_NR

x_NR_m_data_table_4test <- x_NR_m_data_table_R
x_NR_m_data_table_4test$Response <- NULL
names(x_NR_m_data_table_4test) <- c("CNV", "Target_mutation_R", "Others_R")
x_NR_m_data_table_4test$Target_mutation_NR <- x_NR_m_data_table_NR$Target_mutation
x_NR_m_data_table_4test$Others_NR <- x_NR_m_data_table_NR$Others
rownames(x_NR_m_data_table_4test) <- x_NR_m_data_table_4test$CNV 
x_NR_m_data_table_4test$CNV <- NULL
x_NR_m_data_table_4test

# fisher test on each row
row_fisher <- function(row, alt = 'two.sided', cnf = 0.95) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           p_val = f$p.value,
           or = f$estimate[[1]],
           or_ll = f$conf.int[1],
           or_ul = f$conf.int[2]))
}

p <- data.frame(t(apply(x_NR_m_data_table_4test, 1, row_fisher)))
p
write.table(p, "FELINE_gene_cnv_barplot.patient.STable_fisher_test.txt", quote=F, row.name=T, sep="\t")

