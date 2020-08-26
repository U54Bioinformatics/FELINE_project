library(tidyverse)

x= read.table("FEL001046_scRNA.metadata.clinical.txt", sep="\t", header=T)
x$Sample = paste0(x$Sample, "_", x$Timepoint)
x %>% filter(Platform == "10x") %>% select(Cell.ID, Sample, Celltype, Celltype_subtype) -> x_filter
names(x_filter) = c("Cell.ID", "Sample", "Celltype1", "Celltype2")
write.table(x_filter, "FEL011046.cell_type.txt", sep="\t", quote=F, row.names=F, col.names=T)
