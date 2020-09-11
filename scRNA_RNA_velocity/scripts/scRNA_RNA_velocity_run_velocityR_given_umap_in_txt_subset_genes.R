library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
loomfile=args[2]
umap_pos=args[3]
gene_list=args[4]
spliced=args[5]

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

# convert RNA velocity ID back to 10x ID
convert_id_back <- function(ids){
   ids_new = gsub(":", "_", ids)
   ids_new = gsub("x", "", ids_new)
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
umap$Cell.ID = convert_id(umap$Cell.ID)
head(umap)
#row.names(umap) <- umap$Cell.ID
#umap$Cell.ID <- NULL
############################################################################################
##################Load gene list#####################################################
gene_subset=read.table(gene_list, header=T)
gene_subset_id=as.vector(gene_subset$Gene.ID)
######################################################################################


##################RNA velocity analysis######################################################
obj_velocity=paste0(sample, "_velocity.rds")
if (!file.exists(obj_velocity)){
    print("Run new VelocityR analysis")
    # load loom file and convert to seurat obj
    ldat <- ReadVelocity(file =loomfile)
    bm <- as.Seurat(x = ldat)
    # use only filtered cells
    print("umap before overlapping with velocity (bm)")
    print(dim(umap))
    umap <- umap[umap$Cell.ID %in% colnames(bm),]
    print("umap after overlapping with velocity (bm)")
    print(dim(umap))

    print("bm before subseting")
    print(bm)
    #bm <- bm[, colnames(bm) %in% umap$Cell.ID]
    # subset cells
    #bm <- subset(bm, cells = umap$Cell.ID)
    # subset cells and genes
    bm <- subset(bm, cells = umap$Cell.ID, features=gene_subset_id)
    # filter out cells with too few features
    bm <- subset(bm, subset = nFeature_unspliced > 50)

    # use unspliced count
    if (spliced == 'unspliced'){
       DefaultAssay(object = bm) <- "unspliced"
    }
    print("bm after overlap with seurat (pbmc)")
    print(bm)
    # spliced count
    print(bm[["spliced"]])
    print(GetAssayData(object = bm[["spliced"]], slot = "data")[1:5,1:5])
    ## count
    data_to_write_out <- as.data.frame(as.matrix(bm@assays$spliced@data))
    column_name <- names(data_to_write_out)
    names(data_to_write_out) <- convert_id_back(column_name)
    out <- data.table(Gene.ID = rownames(data_to_write_out))
    data_to_write_out <- cbind(out, data_to_write_out)
    fwrite(x = data_to_write_out, row.names = FALSE, sep="\t", file = paste0(sample, '_gene_symbols.spliced.counts.txt'))
    # unspliced count
    print(bm[["unspliced"]])
    print(GetAssayData(object = bm[["unspliced"]], slot = "data")[1:5,1:5])
    data_to_write_out <- as.data.frame(as.matrix(bm@assays$unspliced@data))
    column_name <- names(data_to_write_out)
    names(data_to_write_out) <- convert_id_back(column_name)
    out <- data.table(Gene.ID = rownames(data_to_write_out))
    data_to_write_out <- cbind(out, data_to_write_out)
    fwrite(x = data_to_write_out, row.names = FALSE, sep="\t", file = paste0(sample, '_gene_symbols.unspliced.counts.txt'))
    # ambiguous count
    print(bm[["ambiguous"]])
    print(GetAssayData(object = bm[["ambiguous"]], slot = "data")[1:5,1:5])
    # meta
    out <- data.table(Cell.ID = colnames(x=bm))
    out$Cell.ID <- convert_id_back(out$Cell.ID)
    out_meta <- cbind(out, bm@meta.data)
    fwrite(out_meta, paste(sample, ".cell_metadata.txt", sep=""), sep="\t", quote=F, col.names=T)   
 
    # prepare embeding file
    row.names(umap) <- umap$Cell.ID
    umap$Cell.ID <- NULL
    umap <- as.matrix(umap)

    # perform RNA velocity analysis
    if (spliced == 'spliced'){
       bm <- SCTransform(object = bm, assay = "spliced")
    }
    if (spliced == 'unspliced'){
       bm <- SCTransform(object = bm, assay = "unspliced")
    }
    bm <- RunPCA(object = bm, verbose = FALSE)
    bm <- FindNeighbors(object = bm, dims = 1:20)
    bm <- FindClusters(object = bm)
    bm <- RunUMAP(object = bm, dims = 1:20)
    bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02, spliced.average = 0.05, unspliced.average = 0.05)
    ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
    names(x = ident.colors) <- levels(x = bm)
    cell.colors <- ident.colors[Idents(object = bm)]
    print("color")
    print(head(cell.colors))
    names(x = cell.colors) <- colnames(x = bm)
    print("named color")
    print(head(cell.colors))

    write.table(cell.colors, "cell.colors.txt", sep="\t", row.names=T, quote=F)
    pdf(paste0(sample, "_velocity_denovo_umap.", spliced ,".pdf"))
    #plot umap labeled by clusters
    DimPlot(bm, reduction = "umap", label=TRUE, pt.size = 1)
    #plot umap with RNA velocity
    #default n=200; grid.n=40
    show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm,
       slot = "RunVelocity"), n = 5, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
       cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 4, arrow.lwd = 1,
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
pdf(paste0(sample, "_velocity_given_umap.", spliced ,".pdf"))
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
#ident.colors <- (scales::hue_pal())(n = length(x = levels(x = pbmc)))
#names(x = ident.colors) <- levels(x = pbmc)
#cell.colors <- ident.colors[Idents(object = pbmc)]
#names(x = cell.colors) <- colnames(x = pbmc)

#xcolors = read.table("cell.colors.txt", sep="\t", header=T, comment.char = "")
#cell.colors = as.vector(xcolors$x)
#names(cell.colors) = row.names(xcolors)
#print("color")
#print(head(cell.colors))

print("Plot RNA velocity on given UMAP")
#labeled with cluster, individual cell arrows
vel_out <- show.velocity.on.embedding.cor(emb = umap, vel = Tool(object = bm, slot = "RunVelocity"),
    n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
    cex = 0.8, arrow.scale = 1, show.grid.flow = FALSE, min.grid.cell.mass = 1, grid.n = 40, arrow.lwd = 1,
    do.par = FALSE, cell.border.alpha = 0.1, return.details=F)

#label_clusters(pbmc$seurat_clusters, Embeddings(pbmc, "umap"), font = 2, cex = 1.2)

#labeled with cluster, grid/feild arrows
vel_out <- show.velocity.on.embedding.cor(emb = umap, vel = Tool(object = bm, slot = "RunVelocity"),
    n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
    cex = 0.8, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 1, grid.n = 40, arrow.lwd = 1,
    do.par = FALSE, cell.border.alpha = 0.1, return.details=F)

#label_clusters(pbmc$seurat_clusters, Embeddings(pbmc, "umap"), font = 2, cex = 1.2)


#set color for cells
#ident.colors <- (scales::hue_pal())(n = length(x = levels(x = pbmc)))
#names(x = ident.colors) <- levels(x = pbmc)
#cell.colors <- ident.colors[pbmc$Phase]
#names(x = cell.colors) <- colnames(x = pbmc)
#labeled with cell cycle
#show.velocity.on.embedding.cor(emb = Embeddings(object = pbmc, reduction = "umap"), vel = Tool(object = bm, slot = "RunVelocity"),
#    n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
#    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 1, grid.n = 40, arrow.lwd = 1,
#    do.par = FALSE, cell.border.alpha = 0.1)

#label_clusters(pbmc$Phase, Embeddings(pbmc, "umap"), font = 2, cex = 1.2)

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

print(sessionInfo())

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
