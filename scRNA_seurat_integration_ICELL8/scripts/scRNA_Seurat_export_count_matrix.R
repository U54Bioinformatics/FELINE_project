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
pbmc <- readRDS(paste(prefix, "_Seurat_2kgenes_vst_cc.rds", sep=""))
pbmc <- subset(pbmc, subset = orig.ident != "NA")

#export count for each patient
patients = unique(pbmc$orig.ident)
p_out = data.table(patients)
print(p_out)
#fwrite(x = p_out, row.names = FALSE, col.names=FALSE, sep="\t", file = paste0(prefix, '.patients.list'))
pbmc_small = pbmc

# Write metadata
out <- data.table(Cell.ID = colnames(x=pbmc))
out_meta <- cbind(out, pbmc@meta.data)
fwrite(out_meta, paste(prefix, "cell_metadata.UMAPcluster.txt", sep="."), sep="\t", quote=F, col.names=T)

#for (p in patient) {
#    print(p)
#    pbmc_small = subset(pbmc, subset= orig.ident == p)
    #normalized
    #data_to_write_out <- as.data.frame(as.matrix(pbmc_small@assays$RNA@data))
    #out <- data.table(Gene.ID = rownames(data_to_write_out))
    #data_to_write_out <- cbind(out, data_to_write_out)
    #fwrite(x = data_to_write_out, row.names = FALSE, sep="\t", file = paste0(prefix, '_', p , '_gene_symbols.normalized.counts.txt'))
    #scaled
    #data_to_write_out <- as.data.frame(as.matrix(pbmc_small@assays$RNA@scale.data))
    #out <- data.table(Gene.ID = rownames(data_to_write_out))
    #data_to_write_out <- cbind(out, data_to_write_out)
    #fwrite(x = data_to_write_out, row.names = FALSE, sep="\t", file = paste0(prefix, '_', p , '_gene_symbols.scaled.counts.txt'))
    #raw counts
    mat_raw <- as.data.frame(as.matrix(pbmc_small@assays$RNA@counts))
    id_raw  <- data.table(Gene.ID = rownames(mat_raw))
    #keep raw as int
    mat_raw <- as.data.frame(mat_raw)
    mat_raw <- format(mat_raw, digits=0, nsmall=0)
    out_raw <- cbind(id_raw, mat_raw)
    fwrite(x = out_raw, row.names = FALSE, sep="\t", file = paste0(prefix, '_gene_symbols.raw.counts.txt'))
    #CPM
    mat_cpm <- RelativeCounts(as.matrix(pbmc_small@assays$RNA@counts), scale.factor = 1e6, verbose = TRUE)
    id_cpm <- data.table(Gene.ID = rownames(mat_cpm))
    mat_cpm <- as.data.frame(mat_cpm)
    mat_cpm <- format(mat_cpm, digits=2, nsmall=2)
    out_cpm <- cbind(id_cpm, mat_cpm)
    fwrite(x = out_cpm, row.names = FALSE, sep="\t", file = paste0(prefix, '_gene_symbols.CPM.txt'))
#}#end of for loop
print("Done")
