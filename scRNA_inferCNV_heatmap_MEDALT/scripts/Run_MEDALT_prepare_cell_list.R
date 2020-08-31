library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
patient=args[1]
subclone_dir=args[2]
infercnv_dir=args[3]
medalt_dir=args[4]
expr_rds=args[5]

## filenames
patient_subclone=paste0(subclone_dir, "/", patient, "_HMM_infer.subclone.anno.txt")
patient_cell_final=paste0(medalt_dir, "/", patient, ".cell_final.list")
patient_cell_final_r500=paste0(medalt_dir, "/", patient, ".cell_final.rand500.list")
patient_cnv=paste0(infercnv_dir, "/", patient, ".HMMi6.CNV.txt")
patient_cnv_final=paste0(medalt_dir, "/", patient, ".HMMi6.final.CNV.txt")
patient_cnv_final_names=paste0(medalt_dir, "/", patient, ".HMMi6.final.CNV.names.txt")
patient_cnv_final_r500=paste0(medalt_dir, "/", patient, ".HMMi6.rand500.CNV.txt")
patient_cnv_final_r500_names=paste0(medalt_dir, "/", patient, ".HMMi6.rand500.CNV.names.txt")

## generate cell list for final: which have infercnv results and are retained in final cell list
cell_final <- readRDS(expr_rds)
cell_cnv   <- read.table(patient_subclone, sep="\t", header=T)
cell_cnv_final <- merge(cell_cnv, cell_final,on="Cell.ID")
dim(cell_cnv_final)
cell_cnv_final = data.frame(cell_cnv_final)
cell_cnv_final %>% select(Cell.ID) -> cell_cnv_final
write.table(cell_cnv_final, patient_cell_final, quote=F, row.names=F, col.names=F, sep="\t")
#random 500 cells
try({
write.table(sample(cell_cnv_final$Cell.ID, 500), patient_cell_final_r500, quote=F, row.names=F, col.names=F, sep="\t")})
#sample(as.vector(cell_cnv_final$Cell.ID),100)
#length(cell_cnv_final$Cell.ID)
#cell_cnv_final = cell_cnv_final[,c("Cell.ID")]
#head(cell_cnv_final)

## extract infercnv profile for final cell list
x = read.table(patient_cnv, sep="\t")
dim(x)
x[1:3,1:3]
## final all
y = read.table(patient_cell_final)
x %>% select(as.vector(y$V1)) -> x_out
write.table(names(x_out), patient_cnv_final_names, quote=F, sep="\t", row.names=F, col.names=F)
write.table(x_out, patient_cnv_final, quote=F, sep="\t")
## final rand500
try({
y = read.table(patient_cell_final_r500)
x %>% select(as.vector(y$V1)) -> x_out
write.table(names(x_out), patient_cnv_final_r500_names, quote=F, sep="\t", row.names=F, col.names=F)
write.table(x_out, patient_cnv_final_r500, quote=F, sep="\t")
})

