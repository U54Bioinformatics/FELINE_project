library(Seurat)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
platform=args[2]
#platform="10x"
prefix=paste(sample, platform, sep="_")

#load cell types for UMAP clusters from singleR
pbmc <- readRDS(paste(prefix, "_Seurat_2kgenes_vst_cc.raw.rds", sep=""))
pbmc <- subset(pbmc, subset = orig.ident != "NA")

#apply filters as needed
#pbmc <- subset(pbmc, subset = nFeature_RNA < 7000 & nFeature_RNA > 1000 & nCount_RNA > 2000 & nCount_RNA < 2e5 & percent.mt < 20)

# Write metadata
out <- data.table(Cell.ID = colnames(x=pbmc))
out_meta <- cbind(out, pbmc@meta.data)
fwrite(out_meta, paste(prefix, "cell_metadata.UMAPcluster.txt", sep="_"), sep="\t", quote=F, col.names=T)

# Write count
#mat_raw <- as.data.frame(as.matrix(pbmc_small@assays$RNA@counts))
#id_raw  <- data.table(Gene.ID = rownames(mat_raw))
#keep raw as int
#mat_raw <- as.data.frame(mat_raw)
#mat_raw <- format(mat_raw, digits=0, nsmall=0)
#out_raw <- cbind(id_raw, mat_raw)
#fwrite(x = out_raw, row.names = FALSE, sep="\t", file = paste0(prefix, '_gene_symbols.filtered.counts.txt'))
###Write count
counts <- GetAssayData(object = pbmc, slot = "counts")
counts <- as.matrix(counts)
out <- data.table(Gene.ID = rownames(x=pbmc))
#with Gene.ID
#out_counts <- cbind(out, counts)
#fwrite(out_counts, sep="\t", file = paste0(prefix, '_gene_symbols.filtered.counts.txt'), quote = F)
#without Gene.ID
rownames(counts) <- rownames(x=pbmc)
write.table(counts, row.names = TRUE, sep="\t", file = paste0(prefix, '_gene_symbols.filtered.counts.txt'), quote = F)

print("Done")
