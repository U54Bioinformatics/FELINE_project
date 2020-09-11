library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
loomfile=args[2]
obj_seurat=args[3]

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
       print(id)
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

##################Process preanalyzed Seurat obj############################################
#load pbmc seurat obj
pbmc <- readRDS(obj_seurat)
print("pbmc orignial")
pbmc
head(colnames(pbmc), n = 4)
head(pbmc@meta.data)
#fix cell IDs
#x = gsub(paste0(sample,'_'), paste0(sample,':'), colnames(pbmc))
#x = paste0(x, 'x')
x  = convert_id(colnames(pbmc))
head(x)
length(x)
length(colnames(pbmc))
#colnames(pbmc) <- x
pbmc <- RenameCells(pbmc, new.names=x)
head(colnames(pbmc), n = 4)
############################################################################################

##################RNA velocity analysis######################################################
obj_velocity=paste0(sample, "_velocity.rds")
if (!file.exists(obj_velocity)){
    print("Run new VelocityR analysis")
    #load loom file and convert to seurat obj
    ldat <- ReadVelocity(file =loomfile)
    bm <- as.Seurat(x = ldat)
    print("bm before overlap with pbmc")
    bm
    #use only cell analyzed in seurat obj (filtered cells)
    pbmc <- pbmc[, colnames(pbmc) %in% colnames(bm)]
    print("pbmc after overlap with velocity (bm)")
    pbmc
    bm <- bm[, colnames(bm) %in% colnames(pbmc)]
    print("bm after overlap with seurat (pbmc)")
    bm
    #perform RNA velocity analysis
    bm <- SCTransform(object = bm, assay = "spliced")
    bm <- RunPCA(object = bm, verbose = FALSE)
    bm <- FindNeighbors(object = bm, dims = 1:20)
    bm <- FindClusters(object = bm)
    bm <- RunUMAP(object = bm, dims = 1:20)
    bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
    ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
    names(x = ident.colors) <- levels(x = bm)
    cell.colors <- ident.colors[Idents(object = bm)]
    names(x = cell.colors) <- colnames(x = bm)
    pdf(paste0(sample, "_velocity_denovo_umap.pdf"))
    #plot umap labeled by clusters
    DimPlot(bm, reduction = "umap", label=TRUE, pt.size = 1)
    #plot umap with RNA velocity
    show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm,
       slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
       cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
       do.par = FALSE, cell.border.alpha = 0.1)
    dev.off()
    #save obj
    print("velocity (bm)")
    print(bm)
    head(colnames(bm), n = 4)
    saveRDS(bm, file = obj_velocity)
}else{
    print("Import existing VelocityR analysis")
    bm <- readRDS(obj_velocity)
    print("velocity (bm)")
    print(bm)
    head(colnames(bm), n = 4)
}
############################################################################################


###plot RNA velocity based on UMAP from seurat###########################################################################################################
pdf(paste0(sample, "_velocity_given_umap.pdf"))
#' Get cluster label coordinates
#'
#' @param labels Character or factor vector for labels.
#' @param coords Numeric matrix with two columns, for x and y coordinates of the dimension reduction; the number of rows must match the length of `labels`.
#' @param ... Extra arguments passed to `text`.
#' @return Nothing. Just adds text labels to the plot if the plotting device is still on.
label_clusters <- function(labels, coords, ...) {
  df <- tibble(label = labels, x = coords[,1], y = coords[,2])
  df <- df %>%
    group_by(label) %>%
    summarize(x = median(x), y = median(y))
  text(df$x, df$y, df$label, col="red", ...)
}

#set color for cells as clusters
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = pbmc)))
names(x = ident.colors) <- levels(x = pbmc)
cell.colors <- ident.colors[Idents(object = pbmc)]
names(x = cell.colors) <- colnames(x = pbmc)

#labeled with cluster, individual cell arrows
vel_out <- show.velocity.on.embedding.cor(emb = Embeddings(object = pbmc, reduction = "umap"), vel = Tool(object = bm, slot = "RunVelocity"),
    n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
    cex = 0.8, arrow.scale = 3, show.grid.flow = FALSE, min.grid.cell.mass = 1, grid.n = 40, arrow.lwd = 1,
    do.par = FALSE, cell.border.alpha = 0.1, return.details=F)

label_clusters(pbmc$seurat_clusters, Embeddings(pbmc, "umap"), font = 2, cex = 1.2)

#labeled with cluster, grid/feild arrows
vel_out <- show.velocity.on.embedding.cor(emb = Embeddings(object = pbmc, reduction = "umap"), vel = Tool(object = bm, slot = "RunVelocity"),
    n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 1, grid.n = 40, arrow.lwd = 1,
    do.par = FALSE, cell.border.alpha = 0.1, return.details=F)

label_clusters(pbmc$seurat_clusters, Embeddings(pbmc, "umap"), font = 2, cex = 1.2)


#set color for cells
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = pbmc)))
names(x = ident.colors) <- levels(x = pbmc)
#Idents(object = pbmc) <- 'Phase'
cell.colors <- ident.colors[pbmc$Phase]
names(x = cell.colors) <- colnames(x = pbmc)
#labeled with cell cycle
show.velocity.on.embedding.cor(emb = Embeddings(object = pbmc, reduction = "umap"), vel = Tool(object = bm, slot = "RunVelocity"),
    n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 1, grid.n = 40, arrow.lwd = 1,
    do.par = FALSE, cell.border.alpha = 0.1)

label_clusters(pbmc$Phase, Embeddings(pbmc, "umap"), font = 2, cex = 1.2)

#write.table(vel_out, 'test_velocity.out.txt')

write.table(vel_out$cc, paste0(sample, '_velocity.cc.txt'))
write.table(vel_out$tp, paste0(sample, '_velocity.tp.txt'))
write.table(vel_out$garrows, paste0(sample, '_velocity.garrows.txt'))
write.table(vel_out$arrows, paste0(sample, '_velocity.arrows.txt'))
write.table(vel_out$nd, paste0(sample, '_velocity.nd.txt'))
write.table(vel_out$es, paste0(sample, '_velocity.es.txt'))
write.table(vel_out$gv, paste0(sample, '_velocity.gv.txt'))
write.table(vel_out$gs, paste0(sample, '_velocity.gs.txt'))
write.table(vel_out$scale, paste0(sample, '_velocity.scale.txt'))

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
