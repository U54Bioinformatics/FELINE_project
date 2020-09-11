library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
loomfile=args[2]
umap_pos=args[3]

#loomfile=paste0(fq_dir, '/', sample, '/velocyto/', sample, '.loom')
#run velocity from loom and plot on existing umap
#ldat <- read.loom.matrices(loomfile)
#emat <- ldat$spliced
#nmat <- ldat$unspliced
#read seurat umap coordinate and keep only these cells analyzed by seurat
#seurat_umap <- read.table("VG_CAMA1_10x_Seurat_2kgenes_vst_cc.umap.newid.txt")
#cell.ids  <- rownames(seurat_umap)
#emat <- emat[,cell.ids]
#nmat <- nmat[,cell.ids]
#Estimate RNA velocity (using gene-relative model with k=20 cell kNN pooling and using top/bottom 2% quantiles for gamma fit):
#fit.quantile <- 0.02
#rvel.cd <- gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells=20, fit.quantile=fit.quantile)
#head(rvel.cd)

##################Function to convert cell ID###############################################
#convert 10x ID to velocyto ID
#FEL011_M_1_AATTTCCTCCATCTGC
#FEL011_M_1:AATTTCCTCCATCTGCx
#or
#FEL011_M_AATTTCCTCCATCTGC
#FEL011_M:AATTTCCTCCATCTGCx
#or
#FEL011_AATTTCCTCCATCTGC
#FEL011:AATTTCCTCCATCTGCx
convert_id <- function(ids){
   ids_new = vector()
   for (id in ids){
       #print(id)
       id_array <- strsplit(id, '_')
       index  = id_array[[1]][length(id_array[[1]])]
       string_sample = id_array[[1]][1:length(id_array[[1]])-1]
       sample = paste(string_sample, collapse = '_')
       cell_id_new = paste(sample, index, sep = ":")
       cell_id_new = paste0(cell_id_new, 'x')
       #print(sample)
       #print(index)
       #print(cell_id_new)
       ids_new = c(ids_new, cell_id_new)
   }
   return(ids_new)
}
###########################################################################################

##################Load umap coordinate############################################
#emb <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/embedding.rds"))
#          [,1]       [,2]
#A10  -8.528161  6.5424301
#A11 -17.808458  0.1190365
#A12 -13.631669  3.5772161
umap = read.table(umap_pos, sep="\t", header=T)
#umap$Cell.ID = convert_id(umap$Cell.ID)
head(umap)
#row.names(umap) <- umap$Cell.ID
#umap$Cell.ID <- NULL
#################
#FELINE_test_ArmB.cell_metadata.txt
meta_file = paste0(sample, ".cell_metadata.txt") 
meta = read.table(meta_file, sep="\t", header=T)
merged = merge(umap, meta, by="Cell.ID")
print(head(merged))
############################################################################################
pdf(paste0(sample, "_plot_umap.pdf"))

p = ''
if ( sample %in% c("FEL011046_ArmA", "FELINE_test_ArmA")){
   p <- ggplot(merged, aes(V1, V2))
   xlab <- "V1"
   ylab <- "V2"
}
if ( sample %in% c("FEL011046_ArmB", "FELINE_test_ArmB")){
   p <- ggplot(merged, aes(V1, V3))
   xlab <- "V1"
   ylab <- "V3"
}
if ( sample %in% c("FEL011046_ArmC", "FELINE_test_ArmC")){
   p <- ggplot(merged, aes(V1, V2))
   xlab <- "V1"
   ylab <- "V2"
}

fontsize = 22
p + geom_point(aes(colour=factor(orig.ident))) +
    labs(subtitle="",
       y=ylab,
       x=xlab,
       title="",
       caption = "") +
    theme_classic() + theme_bw()+theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(colour='white'),
                 panel.background = element_blank()) +
   theme(axis.text=element_text(size=fontsize, color='black'),
      #axis.text.x = element_blank(),
      axis.text.x =element_text(size=fontsize, angle = 0, color='black', hjust = 0.8, vjust = 0.8),
      axis.title.x =element_text(size=fontsize, face="bold"),
      axis.title.y =element_text(size=fontsize, face="bold"),
      strip.text = element_text(size=fontsize, face="bold")) +
   theme(legend.text = element_text(size=fontsize),
      legend.title = element_blank(),
      legend.position = "right") +
   theme(plot.title=element_text(size=fontsize-4,face="bold"), axis.text=element_text(size=fontsize, face="bold")) + 
   #geom_vline(xintercept=3*c(1:8)+0.5) +
   #geom_hline(yintercept=1.3) +
   #scale_x_discrete(labels=axis_x_text) +
   #scale_fill_discrete(name = "", labels = unique(x_m_y_z_subset$Subclone1)) +
   theme(plot.margin =  margin(t = 1.5, r = 1, b = 1.5, l = 1, unit = "cm")) +
   guides(colour = guide_legend(override.aes = list(size=8), ncol = 1))

#print(p)
dev.off()

#loomfile="/net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/CellLines/Vince_resistence/preprocess/VG_CAMA1_count03/VG_CAMA1/velocyto/VG_CAMA1.loom"
#ldat <- ReadVelocity(file =loomfile)
#bm <- as.Seurat(x = ldat)
#bm <- SCTransform(object = bm, assay = "spliced")
#bm <- RunPCA(object = bm, verbose = FALSE)
#bm <- FindNeighbors(object = bm, dims = 1:20)
#bm <- FindClusters(object = bm)
#bm <- RunUMAP(object = bm, dims = 1:20)
#bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
#ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
#names(x = ident.colors) <- levels(x = bm)
#cell.colors <- ident.colors[Idents(object = bm)]
#names(x = cell.colors) <- colnames(x = bm)

#pdf(paste0(sample, ".scRNA_velocityR.pdf"))
#plot umap labeled by clusters
#DimPlot(bm, reduction = "umap", label=TRUE, pt.size = 1, group_by = 'ident')
#DimPlot(bm, reduction = "umap", label=TRUE, pt.size = 1, group_by = 'Phase')

#plot umap with RNA velocity
#show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, 
#    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
#    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
#    do.par = FALSE, cell.border.alpha = 0.1)
#dev.off()

#save obj
#saveRDS(bm, file = paste0(sample, ".scRNA_velocityR.rds"))
