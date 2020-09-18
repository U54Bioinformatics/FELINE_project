##this script will read count and meta file from BETSY and create seurat obj for a patient that is specific from agrument
args <- commandArgs(trailingOnly = TRUE)
#sample=args[1]
#platform=args[2]
#patient=args[3]
#prefix=paste(sample, platform, sep="_")
prefix=args[1]
analysis=args[2]
method=args[3]
normalize=as.numeric(args[4])

library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(data.table)
library(RColorBrewer)
library(ggplot2)

featureplot_tsne_png <- function(cells, features, filename, figdir, size=0.5, ncol=3, width=1200, height=1200){
    filename_tsne <- paste(figdir, '/', filename, '.tsne.png', sep="")
    png(filename_tsne, width = width, height = height)
    font_size = 46
    p <- FeaturePlot(cells, features = features, reduction='tsne', slot='data', ncol = ncol, pt.size=size, cols = c("lightgray", "navy", "darkred")) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 28))) + scale_color_gradientn( colours = c('lightgrey', "navy", "darkred"),  limits = c(0, 3))
    print(p)
    dev.off()
}

volinplot_png <- function(cells, features, group, title, filename, figdir, size=2, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 30
    p <- VlnPlot(cells, features = features, group.by = group, pt.size = 0, ncol=3) + labs(title=title, x="", y=features) +  theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8))) + NoLegend() + theme(plot.margin =  margin(t = 1.5, r = 1, b = 1.5, l = 2, unit = "cm"))
    print(p)
    dev.off()
}

dotplot_png <- function(cells, features, group, filename, figdir, size=4, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.dotplot.png', sep="")
    png(filename, width = width, height = height)
    font_size = 46
    p <- DotPlot(cells, features = features, group.by=group) + RotatedAxis() + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    print(p)
    dev.off()
}



pdf(paste(prefix, analysis, "pdf", sep="."))
# Load the PBMC data from 10x count and meta, Lance
counts.df <- fread(paste(prefix, "expression_counts.txt", sep="."), sep="\t", nThread=1)
#counts.df <- as.data.frame(counts.df)
#rownames(counts.df) <- counts.df$Gene.ID
genes <- counts.df$Gene.ID
counts.df <- counts.df[,2:ncol(counts.df)]
counts.df <- as.data.frame(counts.df)
rownames(counts.df) <- genes
print("Count")
dim(counts.df)
counts.df[1:10,1:10]

#write.table(counts, paste(prefix, "gene_symbols.counts.txt", sep="_"), sep="\t", quote = F)
#Annotation
#Annot.df <- read.table(paste0(prefix, ".metadata.txt"), header=TRUE, sep = "\t")
Annot.df <- fread(paste0(prefix, ".metadata.txt"), sep="\t")
print("Annot")
dim(Annot.df)
Annot.df[1:10,1:ncol(Annot.df)]
#count.df <- as.data.frame(counts.df[,Annot.df$Cell.ID])
#rownames(count.df) <- genes
#print("Count cell order")
#dim(count.df)
#count.df
#count_cell <- colnames(count.df)
#head(count_cell)

#rownames(Annot.df) <- Annot.df$Cell.ID
#Annot.df <- Annot.df[,2:ncol(Annot.df)]
#Annot.df  <- Annot.df[match(colnames(counts.df), rownames(Annot.df)),]
##Annot.df <- fread(paste0(prefix, ".metadata.txt"), header=T, sep = "\t", nThread=1)
##Annot.df <- as.data.frame(Annot.df)
##rownames(Annot.df) <- Annot.df$Cell.ID
#count.df <- counts.df[,(colnames(counts.df) %in% Annot.df$Cell.ID)]
##Annot.df <- Annot.df[rownames(Annot.df) %in% colnames(counts.df),]
# reorder Annot.df by counts.df colnames
##count_cell <- colnames(count.df)
##Annot.df <- as.data.frame(Annot.df[count_cell,])
#print("Annot cell order")
#head(Annot.df$Cell.ID)

#print(count.df[1:4,1:4])
#print(Annot.df)
#Create Seurat Object
pbmc <- CreateSeuratObject(counts.df,  project = "T_cells",  meta.data = NULL,  min.cells = 3, min.features = 200)

# Add meta
#celltype <- fread(paste0(prefix, ".metadata.txt"))
rownames(Annot.df) <- Annot.df$Cell.ID
#Annot.df$Cell.ID <- NULL
pbmc <- AddMetaData(pbmc, Annot.df)
###Write metadata
out <- data.table(Cell.ID = colnames(x=pbmc))
out_meta <- cbind(out, pbmc@meta.data)
fwrite(out_meta, paste(prefix, analysis, "cell_metadata.UMAPcluster.txt", sep="."), sep="\t", quote=F, col.names=T)

#Get cell cycle genes
cc.genes <- readLines(con = "/home/jichen/Projects/Breast/scRNA/Data/regev_lab_cell_cycle_genes.txt")
# Segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
#Create scores
#pbmc <- CellCycleScoring(pbmc, s.genes, g2m.genes, set.ident = TRUE)

print(head(pbmc@meta.data))
# add annotation
#celltype <- fread(paste(prefix, ".cell_type.txt", sep=""))
#celltype <- subset(celltype, Cell.ID %in% colnames(pbmc))
#rownames(celltype) <- celltype$Cell.ID
#pbmc <- AddMetaData(pbmc, celltype)

print(head(pbmc@meta.data))

##############################################################################################
###NormalizeData and ScaleData
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(pbmc)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1))
#CombinePlots(plots = list(plot2))
#scaling
#pbmc <- ScaleData(pbmc)
#Regress during scaling
#skip regression first to make cell markers are present then do regression to generate results
if ( TRUE ) {
   print("regression")
   #pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("nCount_RNA", "Sample"))
   #pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("percent.mt", "nCount_RNA","S.Score", "G2M.Score"))
   #pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("percent.mt"))
}
##############################################################################################

#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
#DimPlot(pbmc, reduction = "pca", group.by = "Sample")
#DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#ElbowPlot(pbmc)
if (method == "Seurat"){
   seurat_sds = paste(prefix, analysis, "seurat_standard.rds", sep=".") 
   if (!file.exists(seurat_sds)){
       #seurat
       print("Seurat analysis")
       if (normalize == 1){
           pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
       }
       pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 7000)
       #pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc))
       #pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("Sample"))
       #pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc), vars.to.regress = c("percent.mt", "nCount_RNA", "Sample"))
       #pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"))
       pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nFeature_RNA"))
       pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
       pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:30)
       pbmc <- RunUMAP(pbmc, umap.method = "umap-learn", metric = "correlation", reduction = "pca", dims = 1:30)
       pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:30)
       pbmc <- FindClusters(pbmc, resolution = c(2.0, 1.2, 0.5, 0.8, 1.0))
       saveRDS(pbmc, file = seurat_sds)
   }else{
       print("Existing seurat RDS file: read RDS file")
       pbmc <- readRDS(seurat_sds)
   }
}

if (method == 'FastMNN'){
#FastMNN
print("FastMNN analysis")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- RunFastMNN(object.list = SplitObject(pbmc, split.by = "Sample"))
pbmc <- RunTSNE(pbmc, reduction = "mnn", dims = 1:30)
#pbmc <- RunUMAP(pbmc, umap.method = "umap-learn", metric = "correlation", reduction = "mnn", dims = 1:30)
pbmc <- FindNeighbors(pbmc, reduction = "mnn", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = c(2.0, 1.2, 0.5, 0.8))
}
if (method == "Hormony"){
#Hormony
print("Hormony analysis")
library(harmony)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc, features = rownames(x = pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunHarmony(pbmc, group.by.vars = "Sample")
pbmc <- RunTSNE(pbmc, reduction = "harmony", dims = 1:30) 
pbmc <- RunUMAP(pbmc, umap.method = "umap-learn", metric = "correlation", reduction = "harmony", dims = 1:30)
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = c(1.0, 1.2, 2.0, 0.5, 0.8))
}
if (method == "LIGER"){
#LIGER
print("LIGER analysis")
library(liger)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc, split.by = "Sample", do.center = FALSE)
pbmc <- RunOptimizeALS(pbmc, k = 20, lambda = 5, split.by = "Sample")
pbmc <- RunQuantileAlignSNF(pbmc, split.by = "Sample")
pbmc <- RunTSNE(pbmc, reduction = "iNMF", dims = 1:ncol(pbmc[["iNMF"]]))   
pbmc <- RunUMAP(pbmc, umap.method = "umap-learn", metric = "correlation", dims = 1:ncol(pbmc[["iNMF"]]), reduction = "iNMF")
}

#revised anno
#Annot.df <- fread(paste(prefix, analysis, "cell_metadata.UMAPcluster1.revised_anno.txt", sep="."), sep="\t")
#rownames(Annot.df) <- Annot.df$Cell.ID
#pbmc <- AddMetaData(pbmc, Annot.df)


DimPlot(pbmc, reduction = "tsne", group.by = "RNA_snn_res.0.8", label=TRUE, pt.size = 1) + labs(title = "RNA_snn_res.0.8") 
DimPlot(pbmc, reduction = "tsne", group.by = "RNA_snn_res.1", label=TRUE, pt.size = 1) + labs(title = "RNA_snn_res.1")
DimPlot(pbmc, reduction = "tsne", group.by = "RNA_snn_res.1.2", label=TRUE, pt.size = 1) + labs(title = "RNA_snn_res.1.2")
DimPlot(pbmc, reduction = "tsne", group.by = "RNA_snn_res.2", label=TRUE, pt.size = 1) + labs(title = "RNA_snn_res.2")
DimPlot(pbmc, reduction = "tsne", group.by = "Celltype", label=TRUE, pt.size = 1)
DimPlot(pbmc, reduction = "tsne", group.by = "Celltype_subtype", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "hpca_main_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "Encode_main_type", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "Encode_main_cluster", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "Encode_main_cluster_revised", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "scrublet_doublet_call2", label=TRUE, pt.size = 1)
#DimPlot(pbmc, reduction = "tsne", group.by = "Infercnv_CNA", label=TRUE, pt.size = 1)
DimPlot(pbmc, reduction = "tsne", group.by = "Sample", label=TRUE, pt.size = 1)
DimPlot(pbmc, reduction = "tsne", group.by = "orig.ident", label=TRUE, pt.size = 1)
#Encode_main_type
#Infercnv_CNA

DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.0.8", label=TRUE, pt.size = 1) + labs(title = "RNA_snn_res.0.8")
DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.1", label=TRUE, pt.size = 1) + labs(title = "RNA_snn_res.1")
DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.1.2", label=TRUE, pt.size = 1) + labs(title = "RNA_snn_res.1.2")
DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.2", label=TRUE, pt.size = 1) + labs(title = "RNA_snn_res.2")
DimPlot(pbmc, reduction = "umap", group.by = "Celltype", label=TRUE, pt.size = 1)
DimPlot(pbmc, reduction = "umap", group.by = "Celltype_subtype", label=TRUE, pt.size = 1)
DimPlot(pbmc, reduction = "umap", group.by = "Sample", label=TRUE, pt.size = 1)
DimPlot(pbmc, reduction = "umap", group.by = "orig.ident", label=TRUE, pt.size = 1)
##Write metadata
print("writing metadata")
out <- data.table(Cell.ID = colnames(x=pbmc))
out_meta <- cbind(out, pbmc@meta.data)
fwrite(out_meta, paste(prefix, analysis, "cell_metadata.UMAPcluster1.txt", sep="."), sep="\t", quote=F, col.names=T)

###Write count
if(FALSE){
print("writing countdata")
counts <- GetAssayData(object = pbmc, slot = "counts")
counts <- as.matrix(counts)
out <- data.table(Gene.ID = rownames(x=pbmc))
out_counts <- cbind(out, counts)
fwrite(out_counts, paste(prefix, analysis, "expression_counts.txt", sep="."), sep="\t", quote = F)
}
print("Pipeline finished")
dev.off()

sample=prefix
platform=analysis
figdir=paste(sample, platform, 'png', sep="_")
if(TRUE){
#Immune cell
fig_immune_feature = paste(sample, platform, "Cell_Type_Immunecell_feature", sep="_")
featureplot_tsne_png(pbmc, c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A","CD8B","CD4"), fig_immune_feature, figdir)
fig_immune_violin = paste(sample, platform, "Cell_Type_Immunecell_violin", sep="_")
volinplot_png(pbmc, c("PTPRC", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A","CD8B","CD4"), "seurat_clusters", "Immune", fig_immune_violin, figdir)

#Tcells
fig_tcell_feature = paste(sample, platform, "Cell_Type_Tcell_feature", sep="_")
featureplot_tsne_png(pbmc, c("CD3D", "CD3E", "CD3G", "CD2"), fig_tcell_feature, figdir)
fig_tcell_violin = paste(sample, platform, "Cell_Type_Tcell_violin", sep="_")
volinplot_png(pbmc, c("CD3D", "CD3E", "CD3G", "CD2"), "seurat_clusters", "T cells", fig_tcell_violin, figdir)
#Tcells subset
fig_tcell_feature = paste(sample, platform, "Cell_Type_Tcellsubset_feature", sep="_")
featureplot_tsne_png(pbmc, c("CCR7", "TCF7", "FOXP3", "CD25", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "PD1", "LAG3", "TIGIT", "CTLA4"), fig_tcell_feature, figdir)
fig_tcell_violin = paste(sample, platform, "Cell_Type_Tcellsubset_violin", sep="_")
volinplot_png(pbmc, c("CCR7", "TCF7", "FOXP3", "CD25", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "PD1", "LAG3", "TIGIT", "CTLA4"), "seurat_clusters", "T cells", fig_tcell_violin, figdir)
#Tcells subset2, Analysis of Single-Cell RNA-Seq Identifies Cell-Cell Communication Associated with Tumor Characteristics
fig_tcell_feature = paste(sample, platform, "Cell_Type_Tcellsubset2_feature", sep="_")
featureplot_tsne_png(pbmc, c("CD8A", "CD8B", "CD4", "FOXP3", "IL2RA"), fig_tcell_feature, figdir, ncol=2, width=1000, height=1200)
#fig_tcell_dotplot = paste(sample, platform, "Cell_Type_Tcellsubset2_dotplot_seurat_clusters", sep="_")
#dotplot_png(pbmc, c("CD8A", "CD8B", "CD4", "FOXP3", "IL2RA"), 'seurat_clusters', fig_tcell_dotplot, figdir)
fig_tcell_violin = paste(sample, platform, "Cell_Type_Tcellsubset2_violin", sep="_")
volinplot_png(pbmc, c("CD8A", "CD8B", "CD4", "FOXP3", "IL2RA"), "seurat_clusters", "T cells", fig_tcell_violin, figdir)
#
#fig_tcell_dotplot = paste(sample, platform, "Cell_Type_Tcellsubset2_dotplot_anno", sep="_")
#dotplot_png(pbmc, c("CD8A", "CD8B", "CD4", "FOXP3", "IL2RA"), 'Encode_main_cluster', fig_tcell_dotplot, figdir)
#fig_tcell_violin = paste(sample, platform, "Cell_Type_Tcellsubset2_violin_anno", sep="_")
#volinplot_png(pbmc, c("CD8A", "CD8B", "CD4", "FOXP3", "IL2RA"), 'Encode_main_cluster', "T cells", fig_tcell_violin, figdir)

#PanBcells
fig_panbcell_feature = paste(sample, platform, "Cell_Type_panBcell_feature", sep="_")
featureplot_tsne_png(pbmc, c("CD19", "MS4A1", "CD79A", "CD79B", "CD27", "SDC1"), fig_panbcell_feature, figdir)
fig_panbcell_violin = paste(sample, platform, "Cell_Type_panBcell_violin", sep="_")
volinplot_png(pbmc, c("CD19", "MS4A1", "CD79A", "CD79B", "CD27", "SDC1"), "seurat_clusters", "B cells", fig_panbcell_feature, figdir)
#Bcells
fig_bcell_feature = paste(sample, platform, "Cell_Type_Bcell_feature", sep="_")
featureplot_tsne_png(pbmc, c("CD19", "MS4A1", "CD79A", "CD79B"), fig_bcell_feature, figdir)
fig_bcell_violin = paste(sample, platform, "Cell_Type_Bcell_violin", sep="_")
volinplot_png(pbmc, c("CD19", "MS4A1", "CD79A", "CD79B"), "seurat_clusters", "B cells", fig_bcell_violin, figdir)
#Plasma
fig_plasma_feature = paste(sample, platform, "Cell_Type_plasma_feature", sep="_")
featureplot_tsne_png(pbmc, c("CD27", "SDC1"), fig_plasma_feature, figdir)
fig_plasma_violin = paste(sample, platform, "Cell_Type_plasma_violin", sep="_")
volinplot_png(pbmc, c("CD27", "SDC1"), "seurat_clusters", "Plasma", fig_plasma_violin, figdir)

#Macrophage
#Human monocytes: high levels of CD14 and low CD16 (FCGR3A) and CCR5.
#Human macrophages: low CD14 and high CD16 (FCGR3A), CCR5, CD11b. CD68 and CD71(TFRC) are unique in human macrophages.
#DC: CD1C and CD141(THDB)
#Ben suggestion: CD163, CD206 (MRC1) are M2 marker;CD68 is pan marker;M1=CD68-(CD163+CD206)
fig_macrophage_feature = paste(sample, platform, "Cell_Type_macrophage_feature", sep="_")
featureplot_tsne_png(pbmc, c("CD14", "FCGR3A", "CCR5", "CD68", "TFRC", "CD163", "CSF1R", "CD1C", "MRC1"), fig_macrophage_feature, figdir)
fig_macrophage_violin = paste(sample, platform, "Cell_Type_macrophage_violin", sep="_")
volinplot_png(pbmc, c("CD14", "FCGR3A", "CCR5", "CD68", "TFRC", "CD163", "CSF1R", "CD1C", "MRC1"), "seurat_clusters", "Macrophage", fig_macrophage_violin, figdir)
#M1 and M2
fig_macrophageM12_feature = paste(sample, platform, "Cell_Type_macrophageM12_feature", sep="_")
featureplot_tsne_png(pbmc, c("CD68", "AIF1", "NOS2", "TNF", "IL1B", "CHI3L1", "MRC1", "ARG1", "CD163"), fig_macrophageM12_feature, figdir)
fig_macrophageM12_violin = paste(sample, platform, "Cell_Type_macrophageM12_violin", sep="_")
volinplot_png(pbmc, c("CD68", "AIF1", "NOS2", "TNF", "IL1B", "CHI3L1", "MRC1", "ARG1", "CD163"), "seurat_clusters", "Macrophage", fig_macrophageM12_violin, figdir)
#DC
fig_DC_feature = paste(sample, platform, "Cell_Type_DC_feature", sep="_")
featureplot_tsne_png(pbmc, c("CD14", "CD1C", "THBD", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "ITGAX"), fig_DC_feature, figdir)
fig_DC_violin = paste(sample, platform, "Cell_Type_DC_violin", sep="_")
volinplot_png(pbmc, c("CD14", "CD1C", "THBD", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "ITGAX"), "seurat_clusters", "Macrophage", fig_DC_violin, figdir)

#NK cells
fig_nk_feature = paste(sample, platform, "Cell_Type_NK_feature", sep="_")
featureplot_tsne_png(pbmc, c("NK", "GNLY", "CD56", "CD94", "NKP46", "CD16", "CD57"), fig_nk_feature, figdir)
fig_nk_violin = paste(sample, platform, "Cell_Type_NK_violin", sep="_")
volinplot_png(pbmc, c("NK", "GNLY", "CD56", "CD94", "NKP46", "CD16", "CD57"), "seurat_clusters", "NK cells", fig_nk_violin, figdir)
}

if(TRUE){
#Fibroblast
#Negative marker: PECAM1 (CD31) and KRT1 (cytokeratin)
#Postive marker: ACTA2 (a-SMA), FAP, COL1A1
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_feature", sep="_")
featureplot_tsne_png(pbmc, c("FAP", "COL1A1", "COL3A1", "THY1", "ACTA2", "VIM", "CDH2", "CDH11", "PDPN", "PECAM1", "KRT1"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_violin", sep="_")
volinplot_png(pbmc, c("FAP", "COL1A1", "COL3A1", "THY1", "ACTA2", "VIM", "CDH2", "CDH11", "PDPN", "PECAM1", "KRT1"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
#Fibroblast subset
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_subset_feature", sep="_")
featureplot_tsne_png(pbmc, c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA","FAP", "THY1", "PDPN", "MMP2", "MMP11", "PDGFRL", "TGFB3"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_subset_violin", sep="_")
volinplot_png(pbmc, c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA","FAP", "THY1", "PDPN", "MMP2", "MMP11", "PDGFRL", "TGFB3"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
#Fibroblast subset 2, Fibroblast heterogeneity and immunosuppressive environment in human breast cancer
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_subset2_feature", sep="_")
featureplot_tsne_png(pbmc, c("ITGB1", "FAP", "ACTA2", "PDGFRB", "S100A4", "CAV1", "PECAM1", "KRT1"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_subset2_violin", sep="_")
volinplot_png(pbmc, c("ITGB1", "FAP", "ACTA2", "PDGFRB", "S100A4", "CAV1", "PECAM1", "KRT1"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
#Firboblast secrtome, Cancer-Associated Fibroblasts Their Characteristics and Their Roles in Tumor Growth
#cytokines, chemokines: "ITGB1", "FAP", "ACTA2", "PDGFRB", "S100A4", "CAV1"
#growth factors: EGF, FGF2, IGF1
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_secretome_feature", sep="_")
featureplot_tsne_png(pbmc, c("GDF15", "TGFB2", "CCL5", "CXCL12", "CCL11", "CSF1", "CSF2", "EGF", "FGF2", "IGF1", "IL6"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_secretome_violin", sep="_")
volinplot_png(pbmc, c("GDF15", "TGFB2", "CCL5", "CXCL12", "CCL11", "CSF1", "CSF2", "EGF", "FGF2", "IGF1", "IL6"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
#Fibroblast, collagen, laminin protein in forming extracellular matrix
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_laminin_feature", sep="_")
featureplot_tsne_png(pbmc, c("LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", "COL1A1", "COL1A2"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_laminin_violin", sep="_")
volinplot_png(pbmc, c("LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", "COL1A1", "COL1A2"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)

#epithelial cell
fig_epithelial_feature = paste(sample, platform, "Cell_Type_Epithelial_feature", sep="_")
featureplot_tsne_png(pbmc, c("EPCAM", "ESR1", "ERBB2", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1"), fig_epithelial_feature, figdir)
fig_epithelial_violin = paste(sample, platform, "Cell_Type_Epithelial_violin", sep="_")
volinplot_png(pbmc, c("EPCAM", "ESR1", "ERBB2", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1"), "seurat_clusters", "Epithelial", fig_epithelial_violin, figdir)
}

if(TRUE){
#EMT
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_EMT_feature", sep="_")
featureplot_tsne_png(pbmc, c("CDH1", "CRB3", "DSP", "CDH2", "FN1", "VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_EMT_violin", sep="_")
volinplot_png(pbmc, c("CDH1", "CRB3", "DSP", "CDH2", "FN1", "VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2"), "seurat_clusters", "EMT", fig_fibroblast_violin, figdir)
#NK cells
fig_nk_feature = paste(sample, platform, "Cell_Type_NK_feature", sep="_")
featureplot_tsne_png(pbmc, c("NK", "GNLY", "CD56", "CD94", "NKP46"), fig_nk_feature, figdir)
fig_nk_violin = paste(sample, platform, "Cell_Type_NK_violin", sep="_")
volinplot_png(pbmc, c("NK", "GNLY", "CD56", "CD94", "NKP46"), "seurat_clusters", "NK cells", fig_nk_violin, figdir)

#endothelial
fig_endothelial_feature = paste(sample, platform, "Cell_Type_endothelial_feature", sep="_")
featureplot_tsne_png(pbmc, c("PECAM1", "CD34", "CD36", "ENTPD1", "CD44"), fig_endothelial_feature, figdir)
fig_endothelial_violin = paste(sample, platform, "Cell_Type_endothelial_violin", sep="_")
volinplot_png(pbmc, c("PECAM1", "CD34", "CD36", "ENTPD1", "CD44"), "seurat_clusters", "Endothelial", fig_endothelial_violin, figdir)
#endothelial subset, https://www.nature.com/articles/s41467-018-07770-1
#Pecam1 (CD31), a known pan-endothelial marker, was expressed in both subpopulations.
#Emcn, Cd34 and Sox17, were more abundantly expressed in Vas-Endo
#Lyve1, Prox1, Pdpn, Thy1, and Flt4, were co-expressed in Lym-Endo cells
fig_endothelial_feature = paste(sample, platform, "Cell_Type_endothelial_subset", sep="_")
featureplot_tsne_png(pbmc, c("EMCN", "CD34", "SOX17", "LYVE1", "PROX1", "PDPN", "THY1", "FLT4"), fig_endothelial_feature, figdir)
fig_endothelial_violin = paste(sample, platform, "Cell_Type_endothelial_subset_violin", sep="_")
volinplot_png(pbmc, c("EMCN", "CD34", "SOX17", "LYVE1", "PROX1", "PDPN", "THY1", "FLT4"), "seurat_clusters", "Endothelial", fig_endothelial_violin, figdir)
}
