library(data.table)
library(tidyverse)

# files
file_clinic="FELINE_patient_1_46.clinical_data.txt"
file_subclone="FEL011046_infercnv_subclones_final_HMM_infer.subclone.anno.txt"
file_cellcycle="FEL001046_scRNA.cellcycle.txt"
file_meta="FEL001046_scRNA.metadata.clinical.txt"


# load files
df.clinic = data.frame(fread(file_clinic))
df.clinic = df.clinic[c("Patient.Study.ID", "orig.ident", "Sample", "ARM", "Response")]
df.subclone = data.frame(fread(file_subclone))
names(df.subclone) <- c("Cell.ID", "Subclone")
df.cellcycle = data.frame(fread(file_cellcycle))
df.meta = data.frame(fread(file_meta))
df.meta = df.meta[c("Cell.ID", "Day")]

# merge data
## cell cycle and subclone
df.merged = merge(df.cellcycle, df.subclone, by="Cell.ID")
print("df.merged")
print(dim(df.merged))
print("df.subclone")
print(dim(df.subclone))
print("df.cellcycle")
print(dim(df.cellcycle))

## add clinic
df.merged = merge(df.merged, df.clinic, by="orig.ident")
print("df.merged")
print(dim(df.merged))

## add meta
df.merged = merge(df.merged, df.meta, by="Cell.ID")
print("df.merged")
print(dim(df.merged))
df.merged = df.merged[c("orig.ident", "Patient.Study.ID", "ARM", "Day", "Response", "Cell.ID", "Subclone", "Phase")]

## add genes
gene_list = fread("gene_list.txt", header=FALSE)
gene_data = readRDS("FEL001046_Cancer_cells_scRNA.zinbwave.normalized.RDS")
gene_data %>% filter(Gene.ID %in% gene_list$V1) -> gene_data_sub
row.names(gene_data_sub) <- gene_data_sub$Gene.ID
gene_data_sub$Gene.ID <- NULL
gene_data_sub_t <- t(gene_data_sub)
gene_data_sub_t <- data.frame(gene_data_sub_t)
gene_data_sub_t$Cell.ID <- row.names(gene_data_sub_t)
df.merged_final <- merge(df.merged, gene_data_sub_t, by="Cell.ID")
dim(df.merged_final)


## add pathways
pathway_list = fread("pathway_list.txt", header=FALSE)
pathway_data = readRDS("FEL001046_Cancer_cells_scRNA.zinbwave.normalized.ssGSEA_scores.RDS")
pathway_data %>% filter(`Gene Set` %in% pathway_list$V1) -> pathway_data_sub
row.names(pathway_data_sub) <- pathway_data_sub$`Gene Set`
pathway_data_sub$`Gene Set` <- NULL
pathway_data_sub_t <- t(pathway_data_sub)
pathway_data_sub_t <- data.frame(pathway_data_sub_t)
pathway_data_sub_t$Cell.ID <- row.names(pathway_data_sub_t)
df.merged_final <- merge(df.merged_final, pathway_data_sub_t, by="Cell.ID")
dim(df.merged_final)

## remove genes have NaN, all cells are NaN
sapply(df.merged_final, function(x) sum(is.na(x)))
df.merged_final %>% select_if(~ !any(is.na(.))) -> df.merged_final
fwrite(df.merged_final, "FEL011046_data_gene_pathway.txt", quote=FALSE, sep="\t")


