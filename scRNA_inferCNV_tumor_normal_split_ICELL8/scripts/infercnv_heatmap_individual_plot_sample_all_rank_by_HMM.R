library(ComplexHeatmap)
library(circlize)
library(data.table)
library(fastcluster)
library(RColorBrewer)
set.seed(123)

args <- commandArgs(trailingOnly = TRUE)
#setwd("")
#Create dirs if needed - need to fix
out_loc <- "./"
out_samp <- args[1] #patient name, FEL011
noclust  <- args[2] #0 to skip nocluster, 1 to plot nocluster
k_num    <- args[3] #0 to skip cluster, 2/3../8 to do cluster one each time 
in_type  <- args[4] #infercnv observation file type: pre, hmm, final
#k_num    <- args[1]

print("parameters")
print(out_samp)
print(noclust)
print(k_num)
print(in_type)
print("parameters done")

####Files to load
#Counts
#counts.df <- fread("Expression_Tables/feline_p11.noaggr.counts.txt")
#Chromosome location
chrs <- fread("hg19.RefSeq.NM_pos_unique_sort.txt", header=F)
#Annotation file
annot_df <- fread(paste(out_samp, "_", "HMM_infer.clust", k_num, ".anno.txt", sep = ""), header=T)
names(annot_df) <- c("Cell.ID", "Sample", "nCount_RNA", "Total_Reads_k", "nFeature_RNA", "S", "G2M", "Phase", "percent.mt", "Cell_type", "seurat_clusters", "HMM_inferk")
#colnames(annot_df) <- c("Cell.ID", "Sample")

####Load either CNV or HMM
#Infer CNV Files
infer_txt  = ""
infer_name = ""
if (in_type == "final"){
   infer_txt  = paste(out_samp, ".infercnv.observations.txt", sep="")
   infer_name = "CNA_infer"
} else if (in_type == "hmm"){
   #infer_txt  = "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"
   #infer_txt   = paste(out_samp, ".infercnv.15_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations_dendrogram.observations.txt" , sep="")
   #infer_txt = paste(out_samp, ".infercnv.15_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt", sep="")
   infer_txt = paste(out_samp, ".infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt", sep="")
   #infer_txt  = "infercnv.14_HMM_pred.Bayes_Net.Pnorm_0.5.observations.txt"
   infer_name = "HMM_infer"   
} else if (in_type == "pre"){
   infer_txt  = paste(out_samp, ".infercnv.preliminary.observations.txt", sep="")
   infer_name = "PRE_infer"
}

print(infer_txt)
print(infer_name)

CNA_infer <- read.table(infer_txt, check.name=F, header=T)
out_name  <- infer_name

#CNA_infer <- read.table("infercnv.preliminary.observations.txt", check.name=F, header=T)
#out_name <- "CNA_infer"
#CNA_infer <- read.table("infercnv.12_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations.txt", check.names= F, header=T)
#CNA_infer <- read.table("infercnv.13_HMM_pred.Bayes_Net.Pnorm_0.5.observations.txt", check.names= F, header=T)
#CNA_infer <- read.table("infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt", check.names= F, header=T)
#CNA_infer <- read.table("infercnv.14_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.observations.txt", check.names= F, header=T)
#out_name <- "HMM_infer"
print("check point 1")
#Filter on genes in location file
colnames(chrs) <- c("Gene", "Chr", "Start", "Stop")

#Filter on annotation
#CNA_infer <- CNA_infer[,colnames(CNA_infer) %in% annot_df$Cell.ID]
#annot_df  <- annot_df[match(colnames(CNA_infer), annot_df$Cell.ID),]

#Tranpose CNA matrix
CNA_infer_M <- as.matrix(CNA_infer)
CNA_infer_M <- t(CNA_infer_M)

###Col Annotation
##
#Chromosome top annot InferCNV
chr_sub <- chrs[chrs$Gene %in% colnames(CNA_infer_M), ]
chr_sub$Chr<- factor(chr_sub$Chr, levels= unique(chr_sub$Chr))
chr_Palette <- colorRampPalette(brewer.pal(9, "Set1"))(length(levels(chr_sub$Chr)))
names(chr_Palette) <- levels(chr_sub$Chr)
hc1 = HeatmapAnnotation(Chr = chr_sub$Chr, show_annotation_name=F, show_legend = F,
                        col = list(Chr=chr_Palette))

print("check point 2")
###Row Annotations
##
#Sample
annot_df$Sample<- factor(annot_df$Sample, levels= unique(annot_df$Sample))
#names_Palette <- brewer.pal(length(levels(annot_df$Sample)), "Dark2")
names_Palette <- colorRampPalette(brewer.pal(8, "Dark2"))(length(levels(annot_df$Sample)))
names(names_Palette) <- sort(levels(annot_df$Sample))

hr1 = rowAnnotation(
  Sample = annot_df$Sample, 
  col= list(Sample=names_Palette),
  na_col = "grey", border = TRUE, show_annotation_name = F)

##
#Cell type
annot_df$Cell_type <- factor(annot_df$Cell_type, levels= unique(annot_df$Cell_type))
names_Palette <- colorRampPalette(brewer.pal(8, "Dark2"))(length(levels(annot_df$Cell_type)))
names(names_Palette) <- sort(levels(annot_df$Cell_type))
hr4 = rowAnnotation(
  Encode = annot_df$Cell_type,
  col= list(Encode=names_Palette),
  na_col = "grey", border = TRUE, show_annotation_name = T)

#Seurat UMAP cluster
annot_df$seurat_clusters <- factor(annot_df$seurat_clusters, levels= unique(annot_df$seurat_clusters))
names_Palette <- colorRampPalette(brewer.pal(8, "Spectral"))(length(levels(annot_df$seurat_clusters)))
names(names_Palette) <- sort(levels(annot_df$seurat_clusters))
hr5 = rowAnnotation(
  Cluster = annot_df$seurat_clusters,
  col= list(Cluster=names_Palette),
  na_col = "grey", border = TRUE, show_annotation_name = T)

#Gene Number Hist
hr2 = rowAnnotation(
  Gene_Number = anno_barplot(annot_df$nFeature_RNA, gp = gpar(col = "#666666")),
  na_col = "grey", border = TRUE, show_annotation_name= F)


#Reads Number Hist
#hr3 = rowAnnotation(
#  Read_Number = anno_barplot(annot_df$nCount_RNA, gp = gpar(col = "#666666")),
#  na_col = "grey", border = TRUE, show_annotation_name=F)


#marker genes
counts.df <- fread(paste(out_samp, ".expression_counts.txt", sep=""))
colnames(counts.df)[1] <- "gene_id"
counts.df <- counts.df[counts.df$gene_id %in% chrs$Gene,]
#Genes to add
gene_list <- c("PDPN", "THY1", "CD34", "CDH11", "PTPRC", "GNLY", "CD3E", "CD3D", "CD68", "CD14", "FCGR3A", "CD8A", "ESR1", "EPCAM", "MSLN", "KRT8", "KRT18", "THY1", "CD34", "CDH11")
# Filter counts and log
counts_sub <- counts.df[counts.df$gene_id %in% gene_list, ]
genes <- counts_sub$gene_id
counts_sub <- counts_sub[,colnames(counts_sub) %in% colnames(CNA_infer),with=FALSE]
counts_sub <- setcolorder(counts_sub, colnames(CNA_infer))
counts_sub <- log10(counts_sub +1 )
counts_sub <- as.data.table(t(counts_sub))
colnames(counts_sub) <- genes
#Marker gene colors
col_Fibro = colorRamp2(c(0, 3), c("white", "red"))
col_Imm = colorRamp2(c(0, 3), c("white", "blue"))
col_Epi = colorRamp2(c(0, 3), c("white", "Green"))
#Marker genes
hr3 = rowAnnotation(
  PDPN = as.matrix(counts_sub$PDPN),
  THY1 = as.matrix(counts_sub$THY1),
  CD45 = as.matrix(counts_sub$PTPRC),
  CD3E = as.matrix(counts_sub$CD3E),
  CD68 = as.matrix(counts_sub$CD68),
  EPCAM = as.matrix(counts_sub$EPCAM),
  ESR1 = as.matrix(counts_sub$ESR1),
  KRT8 = as.matrix(counts_sub$KRT8),
  KRT18 = as.matrix(counts_sub$KRT18),
  col = list(PDPN=col_Fibro, THY1=col_Fibro,
             CD45=col_Imm, CD3E=col_Imm, CD68=col_Imm,
             EPCAM=col_Epi, ESR1=col_Epi, KRT8=col_Epi, KRT18=col_Epi),
  border = TRUE,  gap = unit(0, "mm"),
  simple_anno_size = unit(3, "mm"),
  show_legend =rep(FALSE, 9),
  annotation_name_gp =gpar(fontsize=10)
)



# Colors for CNA or
if(out_name == "CNA_infer" || out_name == "PRE_infer"){
  heat_cols = c("navy", "white", "darkred")
} 
if(out_name == "HMM_infer"){
  heat_cols =colorRamp2(c(0, .5, 1, 1.5, 2, 3), c("navy", "royalblue1", "white", "#fee0d2", "#ef3b2c", "#67000d"))
}
print("check point 3")

#order by HMM clustering rank
CNA_infer_M <- CNA_infer_M[annot_df$Cell.ID,]

####No Clustering
if ( as.numeric(noclust) > 0 ) {
    print(c("noclust", noclust))
    print('Plot nocluster plot ......')
    ha <- Heatmap(CNA_infer_M, cluster_rows = F, cluster_columns = F, show_row_names = T, show_heatmap_legend = FALSE,
              show_column_names = F, heatmap_width= unit(20, "cm"), heatmap_height= unit(10, "cm"), col = heat_cols,
              top_annotation = hc1, border=T,
              column_split = factor(chr_sub$Chr, levels= unique(chr_sub$Chr)), cluster_column_slices = FALSE, column_gap = unit(0, "mm"),
              row_split = factor(annot_df$HMM_inferk, levels= sort(unique(annot_df$HMM_inferk))), cluster_row_slices = FALSE, 
              column_title_gp = gpar(fontsize=8), row_title = NULL)

    ht_list <- ha + hr1 + hr2 + hr4 + hr5
    #out for FEL013, too many cell?
    png(paste(out_loc, out_samp, "_", out_name, ".clust_k", k_num ,".png",sep=""), width = 36, height = 12, units = "cm", res=300)
    ht <- draw(ht_list, ht_gap = unit(1, "mm"))
    dev.off()
}#end of if
##############################################block comment out cluster ##########################################


