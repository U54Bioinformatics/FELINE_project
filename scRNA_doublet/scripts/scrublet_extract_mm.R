library(scater)
library(Seurat)
library(scran)
library(BiocSingular)
library(pheatmap)
library(data.table)
library(scds)
library(Matrix)


#the QC process is based on https://f1000research.com/articles/5-2122/v2 and https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/umis.html and https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/doublets.html

args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
timepoint=args[2]
platform=args[3]
#sample="FEL011"
#platform="10x"
prefix=paste(sample, platform, sep="_")

sample_name=paste(sample, timepoint, sep='_')

#read seurat obj and convert to scater obj
pbmc <- readRDS(file = paste0(prefix, "_Seurat_2kgenes_vst_cc.raw.rds"))
pbmc <- subset(pbmc, subset = Sample == sample_name)

###Write count
counts <- GetAssayData(object = pbmc, slot = "counts")
counts <- as.matrix(counts)
#table
write.table(counts, paste(sample_name, ".scrublet.counts.txt", sep=""), sep="\t", quote = F)
#save sparse matrix
sparse.gbm <- Matrix(counts , sparse = T )
head(sparse.gbm)
writeMM(obj = sparse.gbm, file=paste(sample_name, ".matrix.mtx", sep=""))

# save genes and cells names
write(x = rownames(counts), file = paste(sample_name, ".genes.tsv", sep=""))
write(x = colnames(counts), file = paste(sample_name, ".barcodes.tsv", sep=""))

###Write metadata
out <- data.table(Cell.ID = colnames(x=pbmc))
out_meta <- cbind(out, pbmc@meta.data)
fwrite(out_meta, paste(sample_name, ".scrublet.cell_metadata.txt", sep=""), sep="\t", quote=F, col.names=T)

