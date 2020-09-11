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
    p <- FeaturePlot(cells, features = features, reduction='tsne', slot='data', ncol = ncol, pt.size=size, cols = c("lightgray", "navy", "darkred")) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 28))) + scale_color_gradientn( colours = c('lightgrey', "navy", "darkred"))
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

dimplot_png <- function(cells, method, filename, figdir, group, size=4, width=1200, height=1200){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 46
    p <- DimPlot(cells, reduction = method, group.by = group, pt.size=size, label.size=12, label=TRUE) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    print(p)
    dev.off()
}

print_obj_info <- function(pbmc, title){
   print(title)
   print(pbmc)
   try(print(slotNames(pbmc[['RNA']])))
   try(print(slotNames(pbmc[['SCT']])))
   print("counts")
   try(print(GetAssayData(object = pbmc_small[["RNA"]], slot = "counts")[1:5,1:5]))
   print("data")
   try(print(GetAssayData(object = pbmc_small[["RNA"]], slot = "data")[1:5,1:5]))
   print("scale.data")
   try(print(GetAssayData(object = pbmc_small[["RNA"]], slot = "scale.data")[1:5,1:5]))
}


args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
platform=args[2]
umap=0
#sample="FEL011"
#platform="10x"
prefix=paste(sample, platform, sep="_")
figdir=paste(sample, platform, "Seurat_marker_genes", sep="_")

pbmc <- readRDS(paste(prefix, "_Seurat_2kgenes_vst_cc.rds", sep=""))

if(TRUE){
# If the object is generated using SCTransform normalization
# We need to set the assay to RNA
print("Default assay before setting")
print(DefaultAssay(object = pbmc))
DefaultAssay(object = pbmc) <- "RNA"
print("Default assay after setting")
print(DefaultAssay(object = pbmc))
print(slotNames(pbmc[['RNA']]))
# normalized to log1p(CPM)
print("before normalization")
print(slotNames(pbmc[['RNA']]))
try(print(GetAssayData(object = pbmc_small[["RNA"]], slot = "counts")[1:5,1:5]))
try(print(GetAssayData(object = pbmc_small[["RNA"]], slot = "data")[1:5,1:5]))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
print("after normalization")
print(slotNames(pbmc[['RNA']]))
try(print(GetAssayData(object = pbmc_small[["RNA"]], slot = "counts")[1:5,1:5]))
try(print(GetAssayData(object = pbmc_small[["RNA"]], slot = "data")[1:5,1:5]))
}

if(TRUE){
celltype <- fread(paste(prefix, ".broadclass.metadata.txt", sep=""))
rownames(celltype) <- celltype$Cell.ID
pbmc <- AddMetaData(pbmc, celltype)
pbmc <- subset(pbmc, subset = Celltype_subtype != 'Low-quality cells')

fig_tsne = paste(sample, platform, 'tsne', "Celltype", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Celltype")
fig_tsne = paste(sample, platform, 'tsne', "Celltype_subtype", sep="_")
dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Celltype_subtype")
#fig_tsne = paste(sample, platform, 'tsne', "Infercnv_CNA", sep="_")
#dimplot_png(pbmc, 'tsne', fig_tsne, figdir, "Infercnv_CNA")

print('seurat object')
print(pbmc)
print('meta')
print(head(pbmc@meta.data))
}

if(TRUE){
#Single-cell RNA-seq enables comprehensive tumour and immune cell profiling in primary breast Cancer
#Epithelial cells: "KRT19", "CDH1", "EPCAM"
#Stromal cells: "HTRA1", "FBN1", "FAP"
#Immune cells: "PTPRC", "LAPTM5", "IL2RG"
fig_panmarker_feature = paste(sample, platform, "Cell_Type_panmarker_feature", sep="_")
featureplot_tsne_png(pbmc, c("KRT19", "CDH1", "EPCAM", "HTRA1", "FBN1", "FAP", "PTPRC", "LAPTM5", "IL2RG"), fig_panmarker_feature, figdir)
fig_panmarker_violin = paste(sample, platform, "Cell_Type_panmarker_violin", sep="_")
volinplot_png(pbmc, c("KRT19", "CDH1", "EPCAM", "HTRA1", "FBN1", "FAP", "PTPRC", "LAPTM5", "IL2RG"), "seurat_clusters", "Epithelial", fig_panmarker_violin, figdir)
#epithelial cells
fig_epithelial_feature = paste(sample, platform, "Cell_Type_epithelial_feature", sep="_")
featureplot_tsne_png(pbmc, c("KRT19", "CDH1", "EPCAM"), fig_epithelial_feature, figdir)
fig_epithelial_violin = paste(sample, platform, "Cell_Type_epithelial_violin", sep="_")
volinplot_png(pbmc, c("KRT19", "CDH1", "EPCAM"), "seurat_clusters", "Epithelial", fig_epithelial_violin, figdir)
#epithelial cells subset
fig_epithelial_subset_feature = paste(sample, platform, "Cell_Type_epithelial_subset_feature", sep="_")
featureplot_tsne_png(pbmc, c("KRT14", "KRT5", "ACTA2", "MYLK", "TP63", "SLPI", "PROM1", "KRT19", "ANKRD30A", "SYTL2"), fig_epithelial_subset_feature, figdir)
fig_epithelial_subset_violin = paste(sample, platform, "Cell_Type_epithelial_subset_violin", sep="_")
volinplot_png(pbmc, c("KRT14", "KRT5", "ACTA2", "MYLK", "TP63", "SLPI", "PROM1", "KRT19", "ANKRD30A", "SYTL2"), "seurat_clusters", "Epithelial_subset", fig_epithelial_subset_violin, figdir)
#stromal
fig_stromal_feature = paste(sample, platform, "Cell_Type_stromal_feature", sep="_")
featureplot_tsne_png(pbmc, c("HTRA1", "FBN1", "FAP"), fig_stromal_feature, figdir)
fig_stromal_violin = paste(sample, platform, "Cell_Type_stromal_violin", sep="_")
volinplot_png(pbmc, c("HTRA1", "FBN1", "FAP"), "seurat_clusters", "Stromal", fig_stromal_violin, figdir)
#immune
fig_immune_feature = paste(sample, platform, "Cell_Type_immune_feature", sep="_")
featureplot_tsne_png(pbmc, c("PTPRC", "LAPTM5", "IL2RG"), fig_immune_feature, figdir)
fig_immune_violin = paste(sample, platform, "Cell_Type_immune_violin", sep="_")
volinplot_png(pbmc, c("PTPRC", "LAPTM5", "IL2RG"), "seurat_clusters", "Immune", fig_immune_violin, figdir)

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
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_subset3_feature", sep="_")
featureplot_tsne_png(pbmc, c("ACTA2", "TAGLN", "PDGFA", "NMP2", "DCN", "COL1A2", "MYH11"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_subset3_violin", sep="_")
volinplot_png(pbmc, c("ACTA2", "TAGLN", "PDGFA", "NMP2", "DCN", "COL1A2", "MYH11"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
#Cross-species single-cell analysis of pancreatic ductal adenocarcinoma reveals antigen-presenting cancer-associated fibroblasts
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_panCAF_feature", sep="_")
featureplot_tsne_png(pbmc, c("COL1A1", "FAP", "PDPN", "DCN", "VIM"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_panCAF_violin", sep="_")
volinplot_png(pbmc, c("COL1A1", "FAP", "PDPN", "DCN", "VIM"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_myCAF_feature", sep="_")
featureplot_tsne_png(pbmc, c("ACTA2", "TAGLN", "MMP11", "MYL9", "HOPX", "POSTN", "TPM1", "TPM2"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_myCAF_violin", sep="_")
volinplot_png(pbmc, c("ACTA2", "TAGLN", "MMP11", "MYL9", "HOPX", "POSTN", "TPM1", "TPM2"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_iCAF_feature", sep="_")
featureplot_tsne_png(pbmc, c("IL6", "PDGFRA", "CXCL12", "CFD", "DPT", "LMNA", "AGTR1", "HAS1", "CXCL1", "CXCL2", "CCL2", "IL8"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_iCAF_violin", sep="_")
volinplot_png(pbmc, c("IL6", "PDGFRA", "CXCL12", "CFD", "DPT", "LMNA", "AGTR1", "HAS1", "CXCL1", "CXCL2", "CCL2", "IL8"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_apCAF_feature", sep="_")
featureplot_tsne_png(pbmc, c("H2AFB1", "CD74", "SAA3P", "SLPI"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_apCAF_violin", sep="_")
volinplot_png(pbmc, c("H2AFB1", "CD74", "SAA3P", "SLPI"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
#Single cell RNA analysis identifies cellular heterogeneity and adaptive responses of the lung at birth
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_SMC_feature", sep="_")
featureplot_tsne_png(pbmc, c("TGFBI", "PDGFRA", "ACTA2", "ACTG2", "CNN1", "DES"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_SMC_violin", sep="_")
volinplot_png(pbmc, c("TGFBI", "PDGFRA", "ACTA2", "ACTG2", "CNN1", "DES"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_MatrixFB1_feature", sep="_")
featureplot_tsne_png(pbmc, c("TCF21", "WNT2", "FN1", "FGF10", "INMT"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_MatrixFB1_violin", sep="_")
volinplot_png(pbmc, c("TCF21", "WNT2", "FN1", "FGF10", "INMT"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_MatrixFB2_feature", sep="_")
featureplot_tsne_png(pbmc, c("COL1A1", "COL1A2", "AGTR2", "MFAP5", "SFRP2", "IGFBP5"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_MatrixFB2_violin", sep="_")
volinplot_png(pbmc, c("COL1A1", "COL1A2", "AGTR2", "MFAP5", "SFRP2", "IGFBP5"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
#Tissue Specific Origin, Development, and Pathological Perspectives of Pericytes, add CSPG4 (NG2) and MCAM(CD146)
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_pericyte_feature", sep="_")
featureplot_tsne_png(pbmc, c("PDGFRB", "NOTCH3", "MAP3K7CL", "MUSTN1", "ACTA2", "AGTR1A", "VSNL1", "ART3", "CSPG4", "MCAM", "RGS5"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_pericyte_violin", sep="_")
volinplot_png(pbmc, c("PDGFRB", "NOTCH3", "MAP3K7CL", "MUSTN1", "ACTA2", "AGTR1A", "VSNL1", "ART3", "CSPG4", "MCAM", "RGS5"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
#Cancer-associated fibroblast compositions change with breast cancer progression linking S100A4 and PDPN ratios with clinical outcome
#PECAM1 endothelial, RGS5 pericyte, EPCAM epithelial, CD45(PTPRC) immune
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_noCAF_feature", sep="_")
featureplot_tsne_png(pbmc, c("PDPN", "S100A4", "PECAM1", "RGS5", "PECAM", "PTPRC"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_noCAF_violin", sep="_")
volinplot_png(pbmc, c("PDPN", "S100A4", "PECAM1", "RGS5", "PECAM", "PTPRC"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_panCAF1_feature", sep="_")
featureplot_tsne_png(pbmc, c("PDPN", "S100A4"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_panCAF1_violin", sep="_")
volinplot_png(pbmc, c("PDPN", "S100A4"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_pCAF_feature", sep="_")
featureplot_tsne_png(pbmc, c("PDPN", "CXCL12", "FBN1", "SAA3", "ACTA2", "CXCL1", "IL6", "LOX", "COL5A2", "FSTL1", "GPX3"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_pCAF_violin", sep="_")
volinplot_png(pbmc, c("PDPN", "CXCL12", "FBN1", "SAA3", "ACTA2", "CXCL1", "IL6"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_sCAF_feature", sep="_")
featureplot_tsne_png(pbmc, c("HSPD1", "H2AA", "SPP1", "MGP", "CD44", "CLU", "C100A4", "NT5E", "VCAM1"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_sCAF_violin", sep="_")
volinplot_png(pbmc, c("HSPD1", "H2AA", "SPP1", "MGP", "CD44", "CLU", "C100A4", "NT5E", "VCAM1"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
#spatially and functionally distinct subclasses of breast cancer-associated fibroblasts revealed by single cell rna sequencing
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_panCAF2_feature", sep="_")
featureplot_tsne_png(pbmc, c("FAP", "SPARC", "FDGFRA", "FDGFRB", "VIM", "ACTA2"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_panCAF2_violin", sep="_")
volinplot_png(pbmc, c("FAP", "SPARC", "FDGFRA", "FDGFRB", "VIM", "ACTA2"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)
fig_fibroblast_feature = paste(sample, platform, "Cell_Type_fibroblast_pop1_vCAF_feature", sep="_")
featureplot_tsne_png(pbmc, c("ATP1B2", "NOTCH3", "ANO1", "DES", "AOC3", "GUCY1A1", "ESAM", "GDPD3", "MCAM", "HIGD1B", "CPE", "ABCC9", "RGS4", "SPARCL1", "RGS5"), fig_fibroblast_feature, figdir)
fig_fibroblast_violin = paste(sample, platform, "Cell_Type_fibroblast_pop1_vCAF_violin", sep="_")
volinplot_png(pbmc, c("ATP1B2", "NOTCH3", "ANO1", "DES", "AOC3", "GUCY1A1", "ESAM", "GDPD3", "MCAM", "HIGD1B", "CPE", "ABCC9", "RGS4", "SPARCL1", "RGS5"), "seurat_clusters", "Fibroblast", fig_fibroblast_violin, figdir)

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
}

if(TRUE){
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

