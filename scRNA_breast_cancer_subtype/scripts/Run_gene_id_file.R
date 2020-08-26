library(data.table)
library(tidyverse)
counts.df = fread("../scRNA_cellranger/FELINE_cellranger_premRNA/COH049_cell_ranger/read_counts.txt", sep="\t", nThread=1)
counts.df <- subset(counts.df, counts.df$'Gene Symbol' != "")
counts.df <- unique(counts.df, by='Gene Symbol')

counts.df %>% select('Gene ID_1', 'Gene ID_2', 'Gene Symbol') -> x
names(x) = c("Gene.ID_1", "Gene.ID_2", "Gene.Symbol")
write.table(x, "FEL011046.gene_id.txt", quote=F, sep="\t", row.names=F)

