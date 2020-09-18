library(Seurat)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(sctransform)

args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
platform=args[2]
cell.type=args[3]

#fix cell type if these is "_" in it
cell.type.filter = cell.type
if (grepl("_", cell.type)){
    ct_temp = unlist(strsplit(cell.type, '_'))
    cell.type.filter=paste(ct_temp[1], ct_temp[2], sep=" ")
}

prefix=paste(sample, platform, sep="_")

seurat_sds = paste(sample, cell.type, platform, "Seurat_2kgenes_vst_cc.rds", sep="_")
if (!file.exists(seurat_sds)){
    # load cell types for UMAP clusters from singleR
    pbmc <- readRDS(paste(prefix, "_Seurat_2kgenes_vst_cc.rds", sep=""))

    celltype <- fread(paste(prefix, ".metadata.txt", sep=""))
    rownames(celltype) <- celltype$Cell.ID
    pbmc <- AddMetaData(pbmc, celltype)

    # extract cell type obj and meta
    # FEL011043_T_cells_10x_Seurat_2kgenes_vst_cc.rds
    # FEL011043_T_cells_10x.metadata.txt
    pbmc <- subset(pbmc, subset = Celltype == cell.type.filter)
    print("writing rds")
    saveRDS(pbmc, file = paste(sample, cell.type, platform, "Seurat_2kgenes_vst_cc.rds", sep="_"))
}else{
    print("Existing seurat RDS file: read RDS file")
    pbmc <- readRDS(seurat_sds)
}


print("writing metadata")
out <- data.table(Cell.ID = colnames(x=pbmc))
out_meta <- cbind(out, pbmc@meta.data)
fwrite(out_meta, paste(sample, "_", cell.type, "_", platform, ".metadata.txt", sep=""), sep="\t", quote=F, col.names=T)

###Write count
if(TRUE){
print("writing countdata")
DefaultAssay(pbmc) <- 'RNA'
counts <- GetAssayData(object = pbmc, slot = "counts")
counts <- as.matrix(counts)
out <- data.table(Gene.ID = rownames(x=pbmc))
out_counts <- cbind(out, counts)
fwrite(out_counts, paste(sample, "_", cell.type, "_", platform, ".expression_counts.txt", sep=""), sep="\t", quote = F)
}

