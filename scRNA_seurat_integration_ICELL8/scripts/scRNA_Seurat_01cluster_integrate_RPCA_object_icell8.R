args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
platform=args[2]
skip=args[3]
prefix=paste(sample, platform, sep="_")

pdf(paste(prefix, "_Seurat_2kgenes_vst_cc.integrated_RPCA.pdf", sep=""))

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
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA"), group.by="Sample")
VlnPlot(pbmc, features = c("nCount_RNA"), group.by="Sample")
VlnPlot(pbmc, features = c("percent.mt"), group.by="Sample")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="Sample")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Sample")

#filter
#pbmc <- subset(pbmc, subset = nFeature_RNA <= 7000 & nFeature_RNA >= 500 & nCount_RNA >= 500 & nCount_RNA <= 2e7 & percent.mt <= 20)
pbmc <- subset(pbmc, subset = nFeature_RNA < 7000 & nFeature_RNA > 500 & nCount_RNA > 500 & nCount_RNA < 2e7 & percent.mt < 20)
pbmc <- subset(pbmc, subset = Sample != "NA")
#pbmc <- subset(pbmc, subset = orig.ident != "FEL012" & orig.ident != "FEL014" & orig.ident != "FEL015" & orig.ident != "FEL016" & orig.ident != "OB9YER")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA"), group.by="Sample")
VlnPlot(pbmc, features = c("nCount_RNA"), group.by="Sample")
VlnPlot(pbmc, features = c("percent.mt"), group.by="Sample")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="Sample")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Sample")

#Get cell cycle genes
cc.genes <- readLines(con = "/home/jichen/Projects/Breast/scRNA/Data/regev_lab_cell_cycle_genes.txt")
# Segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
#Create scores
pbmc <- CellCycleScoring(pbmc, s.genes, g2m.genes, set.ident = TRUE)

# Remove samples with fewer than 100 cells
x <- table(pbmc$Sample, useNA = c("ifany"))
print("Cell numbers after filter")
x
y <- cbind(names(x), as.numeric(as.character(as.vector(x))))
colnames(y) <- c("Sample", "Cell")
y <- as.data.frame(y)
y$Cell  <- as.numeric(as.character(y$Cell))
y_100 <- subset(y, Cell >= 20)
y_20 <- subset(y, Cell < 20)
y_sample <- as.vector(y_100$Sample)
y_sample20 <- as.vector(y_20$Sample)
print("Sample to keep")
y_sample
#print("Sample to keep: 41 to 43 removed")
#no annotation for 41 to 43, remove 
#y_sample <- y_sample[!y_sample  %in% c("FEL041P137_S", "FEL041P137_M", "FEL041P137_E", "FEL042P401_S", "FEL042P401_M", "FEL042P401_E", "FEL043P202_S", "FEL043P202_M", "FEL043P202_E")]
#y_sample

# output filtered
pbmc_filtered <- pbmc
for (remove_sample in y_sample20){
    pbmc_filtered <- subset(pbmc_filtered, subset = Sample != remove_sample)
}
filtered_obj_file <- paste(prefix, "_Seurat_2kgenes_vst_cc.filtered.rds", sep="")
saveRDS(pbmc_filtered, file = filtered_obj_file)
#meta
out  <- data.table(Cell.ID = colnames(x=pbmc_filtered))
out_meta <- cbind(out,pbmc_filtered@meta.data)
fwrite(out_meta, paste(prefix, "_cell_metadata.UMAPcluster.filtered.txt", sep=""), sep="\t", quote=F, col.names=T)


# Split obj into list and keep only these samples with more than 100 cells
print("Split object")
table(pbmc$Sample)
pbmc.list <- SplitObject(pbmc, split.by = "Sample")
pbmc.list
print("Use subset samples")
pbmc.list <- pbmc.list[y_sample]

#############################AnchorSet#####################################
pbmc.anchors <- ''
AnchorSet_obj_file <- paste(prefix, "_Seurat_2kgenes_vst_cc.integrated_RPCA.AnchorSet.rds", sep="")
if (!file.exists(AnchorSet_obj_file)) {

n_genes = 7000
# Perform Standard normalization separately for each dataset
for (i in 1:length(pbmc.list)) {
    pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], vars.to.regress = c("nCount_RNA"), variable.features.n = n_genes, verbose = TRUE)
    #pbmc.list[[i]] <- NormalizeData(pbmc.list[[i]], verbose = FALSE)
    #pbmc.list[[i]] <- FindVariableFeatures(pbmc.list[[i]], selection.method = "vst", 
    #    nfeatures = 3000, verbose = FALSE)
}

# Next, select features for downstream integration, and run PCA on each object in the list, which is required for running the reciprocal PCA workflow.
#https://github.com/satijalab/seurat/issues/1720
pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = n_genes)
print("pbmc.features")
head(pbmc.features)
print("PrepSCTIntegration")
pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features)
print("RunPCA")
pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
    #x <- ScaleData(x, features = pbmc.features, verbose = FALSE)
    x <- RunPCA(x, features = pbmc.features, verbose = FALSE)
})

# Integrate datasets, and proceed with joint analysis
#https://github.com/satijalab/seurat/issues/997
#This appear to be related to the number of cells in the smallest dataset. 
k.filter <- min(200, min(sapply(pbmc.list, ncol)))
print("FindIntegrationAnchors")
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = pbmc.features, normalization.method = "SCT", reduction = c("rpca"), k.filter = k.filter, reference = c(1, 5, 10), dims = 1:50, verbose = TRUE)
print("pbmc.anchor")
pbmc.anchors
saveRDS(pbmc.anchors, file = AnchorSet_obj_file)
}else{
pbmc.anchors <- readRDS(AnchorSet_obj_file)
}
################################End of AnchorSet##############################################################
Integrated_obj_file = paste(prefix, "_Seurat_2kgenes_vst_cc.integrated_RPCA.rds", sep="")

if (!file.exists(Integrated_obj_file)) {
   pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT", dims = 1:50, verbose = TRUE)
}else{
   pbmc.integrated <- readRDS(Integrated_obj_file)
}
print("integration done")

# add cell type anno
#celltype <- fread(paste(prefix, ".cell_type.txt", sep=""))
#celltype <- subset(celltype, Cell.ID %in% colnames(pbmc))
#rownames(celltype) <- celltype$Cell.ID
#pbmc.integrated <- AddMetaData(pbmc.integrated, celltype)

# Analysis
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(pbmc.integrated) <- "integrated"
#pbmc.integrated <- ScaleData(pbmc.integrated, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE)
pbmc.integrated <- RunPCA(pbmc.integrated, verbose = TRUE)
pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:50, k.param = 20)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = c(2.0, 1.2, 0.5, 0.8))
pbmc.integrated <- RunTSNE(pbmc.integrated, reduction = "pca", dims = 1:50)
pbmc.integrated <- RunUMAP(pbmc.integrated, umap.method = "umap-learn", metric = "correlation", reduction = "pca", dims = 1:50)
DimPlot(pbmc.integrated, group.by = c("Sample"), reduction = "umap", label=TRUE, pt.size = 1, combine = FALSE)
DimPlot(pbmc.integrated, group.by = c("Phase"), reduction = "umap", label=TRUE, pt.size = 1, combine = FALSE)
DimPlot(pbmc.integrated, group.by = c("hpca_main_type"), reduction = "umap", label=TRUE, pt.size = 1, combine = FALSE)
DimPlot(pbmc.integrated, group.by = c("Encode_main_type"), reduction = "umap", label=TRUE, pt.size = 1, combine = FALSE)
DimPlot(pbmc.integrated, group.by = "integrated_snn_res.0.5", reduction = "umap", label=TRUE, pt.size = 1)
DimPlot(pbmc.integrated, group.by = "integrated_snn_res.0.8", reduction = "umap", label=TRUE, pt.size = 1)
DimPlot(pbmc.integrated, group.by = "integrated_snn_res.1.2", reduction = "umap", label=TRUE, pt.size = 1)
DimPlot(pbmc.integrated, group.by = "integrated_snn_res.2", reduction = "umap", label=TRUE, pt.size = 1)
DimPlot(pbmc.integrated, group.by = c("Sample"), reduction = "tsne", label=TRUE, pt.size = 1, combine = FALSE)
DimPlot(pbmc.integrated, group.by = c("Phase"), reduction = "tsne", label=TRUE, pt.size = 1, combine = FALSE)
DimPlot(pbmc.integrated, group.by = c("hpca_main_type"), reduction = "tsne", label=TRUE, pt.size = 1, combine = FALSE)
DimPlot(pbmc.integrated, group.by = c("Encode_main_type"), reduction = "tsne", label=TRUE, pt.size = 1, combine = FALSE)
DimPlot(pbmc.integrated, group.by = "integrated_snn_res.0.5", reduction = "tsne", label=TRUE, pt.size = 1)
DimPlot(pbmc.integrated, group.by = "integrated_snn_res.0.8", reduction = "tsne", label=TRUE, pt.size = 1)
DimPlot(pbmc.integrated, group.by = "integrated_snn_res.1.2", reduction = "tsne", label=TRUE, pt.size = 1)
DimPlot(pbmc.integrated, group.by = "integrated_snn_res.2", reduction = "tsne", label=TRUE, pt.size = 1)

#plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
#CombinePlots(plots)
#p1 + theme(legend.position = "top")
#p2 + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))
#p3 + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))
#p4 + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3)))

dev.off()

###Write metadata
pbmc <-pbmc.integrated
out  <- data.table(Cell.ID = colnames(x=pbmc))
out_meta <- cbind(out, pbmc@meta.data)
fwrite(out_meta, paste(prefix, "_cell_metadata.UMAPcluster.integrated_RPCA.txt", sep=""), sep="\t", quote=F, col.names=T)
#add one columen with cluster# name (numbers are not working with "merge function")
#paste("Cluster", out_meta$seurat_clusters, sep="")
#cluster <- as.data.frame(paste("Cluster", out_meta$seurat_clusters, sep=""))
#colnames(cluster) <- "seurat_clusters2"
#cluster$Cell.ID <- out_meta$Cell.ID
#out_meta <- merge(out_meta, cluster, by="Cell.ID")
#fwrite(out_meta, paste(prefix, "_cell_metadata.UMAPcluster.integrated.txt", sep=""), sep="\t", quote=F, col.names=T)

saveRDS(pbmc, file = paste(prefix, "_Seurat_2kgenes_vst_cc.integrated_RPCA.rds", sep=""))

