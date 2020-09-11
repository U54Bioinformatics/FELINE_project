library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

args <- commandArgs(trailingOnly = TRUE)
fq_dir=args[1]
sample=args[2]
loomfile=paste0(fq_dir, '/veloctyto', sample, '.loom')
#loomfile="/net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/CellLines/Vince_resistence/preprocess/VG_CAMA1_count03/VG_CAMA1/velocyto/VG_CAMA1.loom"
ldat <- ReadVelocity(file =loomfile)
bm <- as.Seurat(x = ldat)
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

pdf(paste0(sample, ".scRNA_velocityR.pdf"))
#plot umap labeled by clusters
DimPlot(bm, reduction = "umap", label=TRUE, pt.size = 1, group_by = 'ident')
DimPlot(bm, reduction = "umap", label=TRUE, pt.size = 1, group_by = 'Phase')

#plot umap with RNA velocity
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, 
    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1)
dev.off()

#save obj
saveRDS(bm, file = paste0(sample, ".scRNA_velocityR.rds"))
