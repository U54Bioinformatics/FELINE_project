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
out_samp <- args[2]
k_num    <- args[1]

####Files to load
#Counts
#counts.df <- fread("Expression_Tables/feline_p11.noaggr.counts.txt")
#Chromosome location
chrs <- fread("/home/jichen/Projects/Breast/scRNA/Data/hg19.RefSeq.NM_pos_unique_sort.txt", header=F)
#Annotation file
annot_df <- fread("cell_metadata.txt", header=T)
names(annot_df) <- c("Cell.ID", "Sample", "nCount_RNA", "Total_Reads_k", "nFeature_RNA", "S", "G2M", "Phase", "percent.mt")
#colnames(annot_df) <- c("Cell.ID", "Sample")

####Load either CNV or HMM
#Infer CNV Files
#CNA_infer <- read.table("infercnv.observations.txt", check.name=F, header=T)
#CNA_infer <- read.table("infercnv.preliminary.observations.txt", check.name=F, header=T)
#out_name <- "CNA_infer"
#CNA_infer <- read.table("infercnv.12_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations.txt", check.names= F, header=T)
#CNA_infer <- read.table("infercnv.13_HMM_pred.Bayes_Net.Pnorm_0.5.observations.txt", check.names= F, header=T)
CNA_infer <- read.table("infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt", check.names= F, header=T)
#CNA_infer <- read.table("infercnv.14_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.observations.txt", check.names= F, header=T)
out_name <- "HMM_infer"

#Filter on genes in location file
colnames(chrs) <- c("Gene", "Chr", "Start", "Stop")

#Filter on annotation
CNA_infer <- CNA_infer[,colnames(CNA_infer) %in% annot_df$Cell.ID]
annot_df  <- annot_df[match(colnames(CNA_infer), annot_df$Cell.ID),]

#Tranpose CNA matrix
CNA_infer_M <- as.matrix(CNA_infer)
CNA_infer_M <- t(CNA_infer_M)

###Col Annotation
##
#Chromosome top annot InferCNV
chr_sub <- chrs[chrs$Gene %in% colnames(CNA_infer_M), ]
chr_sub$Chr<- factor(chr_sub$Chr, levels= unique(chr_sub$Chr))
chr_Palette <- colorRampPalette(brewer.pal(9,"Set1"))(length(levels(chr_sub$Chr)))
names(chr_Palette) <- levels(chr_sub$Chr)
hc1 = HeatmapAnnotation(Chr = chr_sub$Chr, show_annotation_name=F, show_legend = F,
                        col = list(Chr=chr_Palette))


###Row Annotations
##
#Sample
annot_df$Sample<- factor(annot_df$Sample, levels= unique(annot_df$Sample))
names_Palette <- brewer.pal(length(levels(annot_df$Sample)),"Dark2")
names(names_Palette) <- levels(annot_df$Sample)

hr1 = rowAnnotation(
  Sample = annot_df$Sample, 
  col= list(Encode=names_Palette),
  na_col = "grey", border = TRUE, show_annotation_name = F)

#Gene Number Hist
hr2 = rowAnnotation(
  Gene_Number = anno_barplot(annot_df$nFeature_RNA, gp = gpar(col = "#666666")),
  na_col = "grey", border = TRUE, show_annotation_name=F)

#Reads Number Hist
hr3 = rowAnnotation(
  Read_Number = anno_barplot(annot_df$nCount_RNA, gp = gpar(col = "#666666")),
  na_col = "grey", border = TRUE, show_annotation_name=F)

# Colors for CNA or
if(out_name == "CNA_infer"){
  heat_cols = c("navy", "white", "darkred")
} 
if(out_name == "HMM_infer"){
  heat_cols =colorRamp2(c(0, .5, 1, 1.5, 2, 3), c("navy", "royalblue1", "white", "#fee0d2", "#ef3b2c", "#67000d"))
}

####No Clustering
ha <- Heatmap(CNA_infer_M, cluster_rows = F, cluster_columns = F, show_row_names = F, show_heatmap_legend = FALSE,
              show_column_names = F, heatmap_width= unit(20, "cm"), heatmap_height= unit(10, "cm"), col = heat_cols,
              top_annotation = hc1, border=T,
              column_split = factor(chr_sub$Chr, levels= unique(chr_sub$Chr)), cluster_column_slices = FALSE, column_gap = unit(0, "mm"),
              row_split = factor(annot_df$Sample, levels= unique(annot_df$Sample)), cluster_row_slices = FALSE, 
              column_title_gp = gpar(fontsize=8), row_title = NULL)

ht_list <- ha + hr1 + hr2 + hr3

#png(paste(out_loc, out_samp, out_name, "_NoClust.png",sep=""), width = 36, height = 12, units = "cm", res=300)
#ht <- draw(ht_list, ht_gap = unit(1, "mm"))
#dev.off()

#Cluster method
dend <- fastcluster::hclust(dist(CNA_infer_M), method="ward.D2")
###############for loop for k#######################################################################
for (k in as.numeric(k_num)){
####K5
ha <- Heatmap(CNA_infer_M, cluster_rows = dend, cluster_columns = F, show_row_names = F, show_heatmap_legend = FALSE,
              show_column_names = F, heatmap_width= unit(20, "cm"), heatmap_height= unit(10, "cm"), col = heat_cols,
              clustering_method_rows = "ward.D2", top_annotation = hc1,
              row_split = k, row_gap = unit(0, "mm"), border=T,
              column_split = factor(chr_sub$Chr, levels= unique(chr_sub$Chr)), cluster_column_slices = FALSE, column_gap = unit(0, "mm"),
              column_title_gp = gpar(fontsize=8))

ht_list <- ha + hr1 + hr2 + hr3
png(paste(out_loc, out_samp, out_name, "_k", k, ".png",sep=""), width = 36, height = 12, units = "cm", res=300)
ht = draw(ht_list, ht_gap = unit(1, "mm"))
dev.off()

#Get clusters
out_list <- row_order(ht)
names(out_list) <- 1:k

for (i in 1:length(out_list)){
  if (i == 1) {
    print(i)
    print(names(out_list)[i][1])
    temp <- row.names(CNA_infer_M[out_list[[i]],])
    out <- cbind(temp, paste("Cluster", names(out_list)[i][1], sep=""))
    colnames(out) <- c("Cell.ID", "Cluster")
  } 
  else {
    print(i)
    print(names(out_list)[i][1])
    clu <- row.names(CNA_infer_M[out_list[[i]],])
    clu <- cbind(clu, paste("Cluster", names(out_list)[i][1], sep=""))
    colnames(clu) <- c("Cell.ID", "Cluster")
    out <- rbind(out, clu)
  }
}
#Annot k5
colnames(out) <- c("Cell.ID", paste(out_name, "k", k, sep=""))
out <- as.data.table(out)
annot_df <- merge(annot_df, out, by = "Cell.ID")
}###############for loop for k#####################################################################

####K4##########block comment######################################################################
if (FALSE){
ha <- Heatmap(CNA_infer_M, cluster_rows = dend, cluster_columns = F, show_row_names = F, show_heatmap_legend = FALSE,
              show_column_names = F, heatmap_width= unit(20, "cm"), heatmap_height= unit(10, "cm"), col = heat_cols,
              clustering_method_rows = "ward.D2", top_annotation = hc1,
              row_split = 8, row_gap = unit(0, "mm"), border=T,
              column_split = factor(chr_sub$Chr, levels= unique(chr_sub$Chr)), cluster_column_slices = FALSE, column_gap = unit(0, "mm"),
              column_title_gp = gpar(fontsize=8))

ht_list <- ha + hr1 
png(paste(out_loc, out_samp, out_name, "_k8.png",sep=""), width = 36, height = 12, units = "cm", res=300)
ht = draw(ht_list, ht_gap = unit(1, "mm"))
dev.off()

out_list <- row_order(ht)
names(out_list) <- c("1", "2", "3", "4", "5", "6", "7", "8")

#Get clusters
for (i in 1:length(out_list)){
  if (i == 1) {
    print(i)
    print(names(out_list)[i][1])
    temp <- row.names(CNA_infer_M[out_list[[i]],])
    out <- cbind(temp, paste("cluster", names(out_list)[i][1], sep=""))
    colnames(out) <- c("Cell.ID", "Cluster")
  } 
  else {
    print(i)
    print(names(out_list)[i][1])
    clu <- row.names(CNA_infer_M[out_list[[i]],])
    clu <- cbind(clu, paste("cluster", names(out_list)[i][1], sep=""))
    colnames(clu) <- c("Cell.ID", "Cluster")
    out <- rbind(out, clu)
  }
}
#Annot k4
colnames(out) <- c("Cell.ID", paste(out_name, "k8",sep=""))
out <- as.data.table(out)
annot_df <- merge(annot_df, out, by = "Cell.ID")
}################block comment##################################################################

fwrite(annot_df, paste(out_samp, out_name, ".clust", k_num, ".anno.txt", sep = ""), sep="\t", quote = F)




