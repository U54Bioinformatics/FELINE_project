library(dplyr)
library(Seurat)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
sample=args[1]
platform=args[2]
control_count=args[3]
control_annot=args[4]
#skip=args[3]
#sample="FEL011"
#platform="10x"
prefix=paste(sample, platform, sep="_")

#Load exp and wellList
counts.df <-  fread(paste(sample, platform, "count01.noaggr.counts.txt", sep="_"), sep="\t", nThread=1)
#counts.df <- subset(counts.df, counts.df$'Gene Symbol' != "")
#counts.df <- unique(counts.df, by='Gene Symbol')
counts.df <- as.data.frame(counts.df)
rownames(counts.df) <- counts.df$`gene_id`
counts.df <- counts.df[,3:ncol(counts.df)]

#write.table(counts, paste(prefix, "gene_symbols.counts.txt", sep="_"), sep="\t", quote = F)
#Annotation
Annot.df <- read.table(paste(sample, platform, "cell_metadata.txt", sep="_"), header=TRUE, sep = "\t")
rownames(Annot.df) <- Annot.df$Cell
#Annot.df <- Annot.df[,2:ncol(Annot.df)]
Annot.df <- Annot.df[,1:2]
Annot.df <- Annot.df[match(colnames(counts.df), rownames(Annot.df)),]

#read cell line matrix
Fibro.data <- read.table(control_count, check.names=FALSE)
colnames(Fibro.data) <- colnames(Fibro.data)
dim(Fibro.data)
cellline.data <- Fibro.data

#read cell line anno
cellline.anno <- fread(control_annot, sep="\t", nThread=1)
rownames(cellline.anno) <- cellline.anno$Cell
anno_merged <- rbind(Annot.df, cellline.anno)
write.table(anno_merged, file=paste(prefix, "_counts.anno.raw.txt", sep=""), quote=F, row.names=FALSE, col.names=FALSE, sep="\t")
#cellline.anno <- cellline.anno[,2:ncol(cellline.anno)]
#anno_merged <- merge(Annot.df, cellline.anno, by = "Sample", all = TRUE)
#anno_merged$Cell.ID <- NULL

#merged cell line data
#cellline.data <- merge(Fibro.data, HMEC.data, by = "row.names", all = TRUE)
#rownames(cellline.data) <- cellline.data$Row.names
#cellline.data <- cellline.data[-1]
#dim(cellline.data)

#merge these two datasets
pbmc_merged <- merge(counts.df, cellline.data, by = "row.names", all = TRUE)
#pbmc_merged[is.na(pbmc_merged)] <- 0
row.has.na <- apply(pbmc_merged, 1, function(x){any(is.na(x))})
pbmc_merged <- pbmc_merged[!row.has.na,]
rownames(pbmc_merged) <- pbmc_merged$Row.names
pbmc_merged <- pbmc_merged[-1]
dim(pbmc_merged)

#write merged data to file
write.table(pbmc_merged, file=paste(prefix, "_counts.raw.matrix", sep=""), quote=F, sep="\t", row.names=TRUE)


if ( FALSE ) {
pdf(paste(prefix, "_Seurat_2kgenes_vst_cc.pdf", sep=""))

#annotation
pbmc <- CreateSeuratObject(counts = pbmc_merged, meta=anno_merged, project = sample, min.cells = 3, min.features = 200)
#x <- paste(pbmc$Cell, pbmc$Sample, sep="\t")
#write.table(x, file=paste(prefix, "_counts.anno.raw.txt", sep=""), quote=F, col.names=FALSE, sep="\t")

dim(counts.df)
dim(cellline.data)
dim(pbmc_merged)

#
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA"), ncol = 1)
VlnPlot(pbmc, features = c("nCount_RNA"), ncol = 1)
VlnPlot(pbmc, features = c("percent.mt"), ncol = 1)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 50)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))

#scaling
pbmc <- ScaleData(pbmc)
#Get cell cycle genes
cc.genes <- readLines(con = "/home/jichen/Projects/Breast/scRNA/Data/regev_lab_cell_cycle_genes.txt")
# Segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
#Create scores
pbmc <- CellCycleScoring(pbmc, s.genes, g2m.genes, set.ident = TRUE)
#Regress during scaling
#pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("percent.mt", "nCount_RNA","S.Score", "G2M.Score"))
#pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("percent.mt"))

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:25)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:25)

#plot features on top of cluster and cells
FeaturePlot(object = pbmc, features = c("nCount_RNA"), pt.size = 4)
FeaturePlot(object = pbmc, features = c("nFeature_RNA"), pt.size = 4)

DimPlot(pbmc, reduction = "umap", label=TRUE)
#labeled clustered by samples (orig.ident)
DimPlot(pbmc, reduction = "umap", group.by = "orig.ident", label=TRUE)

#Immune cell
FeaturePlot(pbmc, features = c("PTPRC", "CD68", "CD14", "FCGR3A"))
VlnPlot(pbmc, features = c("PTPRC", "CD68", "CD14", "FCGR3A"))
#Fibroblast cell
FeaturePlot(pbmc, features = c("PDPN", "CD34"))
VlnPlot(pbmc, features = c("PDPN", "CD34"))
#epithelial cell
FeaturePlot(pbmc, features = c("EPCAM", "ESR1", "KRT7", "KRT8", "KRT18", "KRT19"))
VlnPlot(pbmc, features = c("EPCAM", "ESR1", "KRT7", "KRT8", "KRT18", "KRT19"))

dev.off()
}



