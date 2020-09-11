#export R_LIBS=/home/jichen/software/BETSY/install/envs/zinbwave/lib/R/library
library(umap)
library(data.table)
library(dplyr)

#ARM A full data
infile="Jason_input/10_umap_visualization/ARM_A_cancer_cell_UMAP_results.RData"
load(file=infile)
# umap position
umap_cor = u_datFULL[,c("Cell.ID", "V1", "V2")]
row.names(umap_cor) <- umap_cor$Cell.ID
fwrite(umap_cor, file="Jason_input/10_umap_visualization/ARM_A_cancer_cell_UMAP_results.umap.txt", sep="\t", quote=FALSE, row.names=FALSE)
# gene list
gene_list_a = data.table(ss_dd%>%colnames)
names(gene_list_a) <- "Gene.ID"
fwrite(gene_list_a, file="Jason_input/10_umap_visualization/ARM_A_cancer_cell_UMAP_results.gene.txt", sep="\t", quote=FALSE, row.names=FALSE)

#ARM B full data
infile="Jason_input/10_umap_visualization/ARM_B_Day180_cancer_cell_UMAP_results.RData"
load(file=infile)
umap_cor = u_datFULL[,c("Cell.ID", "V1", "V3")]
row.names(umap_cor) <- umap_cor$Cell.ID
fwrite(umap_cor, file="Jason_input/10_umap_visualization/ARM_B_Day180_cancer_cell_UMAP_results.umap.txt", sep="\t", quote=FALSE, row.names=FALSE)
# gene list
gene_list_b = data.table(ss_dd%>%colnames)
names(gene_list_b) <- "Gene.ID"
fwrite(gene_list_b, file="Jason_input/10_umap_visualization/ARM_B_Day180_cancer_cell_UMAP_results.gene.txt", sep="\t", quote=FALSE, row.names=FALSE)

#ARM C full data
infile="Jason_input/10_umap_visualization/ARM_C_Day180_cancer_cell_UMAP_results.RData"
load(file=infile)
umap_cor = u_datFULL[,c("Cell.ID", "V1", "V2")]
row.names(umap_cor) <- umap_cor$Cell.ID
fwrite(umap_cor, file="Jason_input/10_umap_visualization/ARM_C_Day180_cancer_cell_UMAP_results.umap.txt", sep="\t", quote=FALSE, row.names=FALSE)
# gene list
gene_list_c = data.table(ss_dd%>%colnames)
names(gene_list_c) <- "Gene.ID"
fwrite(gene_list_c, file="Jason_input/10_umap_visualization/ARM_C_Day180_cancer_cell_UMAP_results.gene.txt", sep="\t", quote=FALSE, row.names=FALSE)

