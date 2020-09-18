args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
platform=args[2]
skip=args[3]
prefix=paste(sample, platform, sep="_")

library(dplyr)
library(Seurat)
library(data.table)
library(ggplot2)
#do not use future and "multiprocess", will stunk
# set up future for parallelization
#library(future)
#library(future.apply)
#plan("multiprocess", workers = 4)
#options(future.globals.maxSize = 100000 * 1024^2)
# Load the PBMC dataset from datadir
#pbmc.data <- Read10X(data.dir = paste(prefix, "_datadir/filtered_gene_bc_matrices/hg19", sep=""), gene.column=2)
# Initialize the Seurat object with the raw (non-normalized data).
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = sample, min.cells = 3, min.features = 200)
# merge obj into a single obj
pbmc <- ''
merged_obj_file <- paste(prefix, "_Seurat_2kgenes_vst_cc.merged.rds", sep="")
if (!file.exists(merged_obj_file)) {
   FEL001_obj <- readRDS('FEL001P120_icell8_Seurat_2kgenes_vst_cc.raw.rds')
   FEL002_obj <- readRDS('FEL002P106_icell8_Seurat_2kgenes_vst_cc.raw.rds')
   FEL003_obj <- readRDS('FEL003P112_icell8_Seurat_2kgenes_vst_cc.raw.rds')
   FEL004_obj <- readRDS('FEL004P114_icell8_Seurat_2kgenes_vst_cc.raw.rds')
   FEL005_obj <- readRDS('FEL005P133_icell8_Seurat_2kgenes_vst_cc.raw.rds')
   FEL006_obj <- readRDS('FEL006P111_icell8_Seurat_2kgenes_vst_cc.raw.rds')
   FEL007_obj <- readRDS('FEL007P127_icell8_Seurat_2kgenes_vst_cc.raw.rds')
   FEL008_obj <- readRDS('FEL008P108_icell8_Seurat_2kgenes_vst_cc.raw.rds')
   FEL009_obj <- readRDS('FEL009P136_icell8_Seurat_2kgenes_vst_cc.raw.rds')
   FEL010_obj <- readRDS('FEL010P104_icell8_Seurat_2kgenes_vst_cc.raw.rds')
   pbmc <- merge(FEL001_obj, y = c(FEL002_obj, FEL003_obj, FEL004_obj, FEL005_obj, FEL006_obj, FEL007_obj, FEL008_obj, FEL009_obj, FEL010_obj), add.cell.ids = c("FEL001", "FEL002", "FEL003", "FEL004", "FEL005", "FEL006", "FEL007", "FEL008", "FEL009", "FEL010"), project = sample)
   #meta
   out  <- data.table(Cell.ID = colnames(x=pbmc))
   out_meta <- cbind(out, pbmc@meta.data)
   fwrite(out_meta, paste(prefix, "_cell_metadata.UMAPcluster.merged.txt", sep=""), sep="\t", quote=F, col.names=T)
   #obj
   saveRDS(pbmc, file = merged_obj_file)
}else{
   pbmc <- readRDS(merged_obj_file) 
}

#pbmc <- subset(pbmc, subset = orig.ident != "FEL041" & orig.ident != "FEL042" & orig.ident != "FEL043")

head(colnames(pbmc))
tail(colnames(pbmc))
unique(sapply(X = strsplit(colnames(pbmc), split = "_"), FUN = "[", 1))
table(pbmc$orig.ident, useNA = c("ifany"))
print("Cell numbers before filter")
table(pbmc$Sample, useNA = c("ifany"))

# add annotation
celltype <- fread(paste(prefix, ".cell_type.txt", sep=""))
celltype <- subset(celltype, Cell.ID %in% colnames(pbmc))
rownames(celltype) <- celltype$Cell.ID
pbmc <- AddMetaData(pbmc, celltype)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#filter
#pbmc <- subset(pbmc, subset = nFeature_RNA <= 7000 & nFeature_RNA >= 500 & nCount_RNA >= 500 & nCount_RNA <= 2e7 & percent.mt <= 20)
pbmc <- subset(pbmc, subset = nFeature_RNA < 7000 & nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 2e7 & percent.mt < 20)
pbmc <- subset(pbmc, subset = Sample != "NA")
pbmc <- subset(pbmc, subset = orig.ident != "")
#pbmc <- subset(pbmc, subset = orig.ident != "FEL012" & orig.ident != "FEL014" & orig.ident != "FEL015" & orig.ident != "FEL016" & orig.ident != "OB9YER")

pbmc_filtered <- pbmc
filtered_obj_file <- paste(prefix, "_Seurat_2kgenes_vst_cc.filtered.rds", sep="")
saveRDS(pbmc_filtered, file = filtered_obj_file)
#meta
out  <- data.table(Cell.ID = colnames(x=pbmc_filtered))
out_meta <- cbind(out,pbmc_filtered@meta.data)
fwrite(out_meta, paste(prefix, "_cell_metadata.UMAPcluster.filtered.txt", sep=""), sep="\t", quote=F, col.names=T)

