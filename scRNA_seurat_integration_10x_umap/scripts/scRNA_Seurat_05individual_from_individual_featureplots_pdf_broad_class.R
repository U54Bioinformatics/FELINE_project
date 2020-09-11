library(Seurat)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

dimplot_png <- function(cells, method, filename, figdir, group, size=4, width=1200, height=1000){
    filename <- paste(figdir, '/', filename, '.png', sep="")
    png(filename, width = width, height = height)
    font_size = 46
    p <- DimPlot(cells, reduction = method, group.by = group, pt.size=size, label.size=12, label=TRUE) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
    print(p)
    dev.off()
}

dotplot_png <- function(cells, features, group, filename, figdir, size=4, width=1200, height=1200, pdf){
    if (as.numeric(pdf)){
        filename <- paste(figdir, '/', filename, '.dotplot_pdf.pdf', sep="")
        pdf(filename, width = width, height = height)
        font_size = 26
        p <- DotPlot(cells, features = features, group.by=group, dot.scale = 16) + RotatedAxis() + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
        print(p)
        dev.off()
    }else{
        filename <- paste(figdir, '/', filename, '.dotplot_png.png', sep="")
        png(filename, width = width, height = height)
        font_size = 46
        p <- DotPlot(cells, features = features, group.by=group, dot.scale = 16) + RotatedAxis() + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + guides(shape = guide_legend(override.aes = list(size = 8)))
        print(p)
        dev.off()
    }
}

args <- commandArgs(trailingOnly = TRUE)
sample=args[1]
platform=args[2]
umap=0
#sample="FEL011"
#platform="10x"
prefix=paste(sample, platform, sep="_")
figdir=paste(sample, platform, "Seurat_UMAP_Celltype", sep="_")

if (TRUE){
    pbmc <- readRDS(paste(prefix, "_Seurat_2kgenes_vst_cc.rds", sep=""))
    DefaultAssay(object = pbmc) <- "RNA"

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

        # dotplot for marker genes
        # Epithelial: 'KRT19', 'CDH1'; Fibroblasts: 'FAP', 'COL1A1'; Endothelial: 'PECAM1'; Pericytes: 'PDGFRB', 'NOTCH3'
        # Immuune: 'PTPRC'; Macrophages: 'CSF1R'; T cells: 'CD4', 'CD8A', 'CD8B', B cells: 'MS4A1'; Plasma cells: 'SDC1'
        fig_tsne = paste(sample, platform, 'tsne', "marker_genes_dotplot", sep="_")
        markers  = c('KRT19', 'CDH1', 'FAP', 'COL1A1', 'PECAM1', 'PDGFRB', 'NOTCH3', 'PTPRC', 'CSF1R', 'CD3D', 'MS4A1', 'SDC1')
        dotplot_png(pbmc, markers, 'Celltype', fig_tsne, figdir, size=4, width=18, height=8, 1)
        dotplot_png(pbmc, markers, 'Celltype', fig_tsne, figdir, size=4, width=1800, height=800, 0)
        print('seurat object')
        print(pbmc)
        print('meta')
        print(head(pbmc@meta.data))
    }

    revised_head <- 'Celltype'
    fig_tsne = paste(sample, platform, revised_head, sep="_")

    print('Extract embedding')
    coords <- as.data.frame(Embeddings(object =pbmc, reduction = 'tsne'))
    rownames(celltype) <- celltype$Cell.ID
    pbmc <- AddMetaData(pbmc, celltype)
    Idents(pbmc) <- pbmc[[revised_head]]

    print('coords and celltype')
    head(coords)
    head(celltype)

    print('TSNE_1 and TSNE_2')
    pbmc$tSNE_1 <- coords$tSNE_1
    pbmc$tSNE_2 <- coords$tSNE_2
    #x <- cbind(pbmc$UMAP_1, pbmc$UMAP_2, pbmc$Encode_main_cluster)
    #x <- cbind(pbmc$UMAP_1, pbmc$UMAP_2, pbmc$seurat_clusters)
    x <- cbind(pbmc$tSNE_1,pbmc$tSNE_2,Idents(pbmc))
    colnames(x) <- c("tSNE_1", "tSNE_2", "Cell.Type")
    x <- as.data.frame(x)

    # output TSNE coordinate
    #x1 <- cbind(pbmc$Cell.ID, pbmc$tSNE_1,pbmc$tSNE_2,Idents(pbmc))
    #colnames(x) <- c("Cell.ID", "tSNE_1", "tSNE_2", "Cell.Type")
    #fwrite(x, paste(sample, platform, "TSNE_coordinate.txt", sep="_"), sep="\t", quote=F, col.names=T)
}else{
    print('pass')
    #x1 <- fread(paste(sample, platform, "TSNE_coordinate.txt", sep="_"), sep="\t")
    #x1$Cell.ID  <- NULL
    #x  <- x1
}

    #cell type -> color
    c <- data.frame(
            Cell.Type = c("Cancer cells", "Normal epithelial cells", "Stromal cells", "Immune cells"),
            Color = c("#E7298A", "#66A61E", "#006633", "light sea green")
    )
    c1 <- data.frame(
            Cell.Type = c("Cancer cells", "Adipocytes", "B cells", "T cells", "Endothelial cells", "Normal epithelial cells", "Fibroblasts", "Macrophages", "Plasma cells", "Pericytes", "Low-quality cells"),
            Color = c("#E7298A", "#7570B3", "#994C00", "#CC99FF", "#FF9933", "#66A61E", "#006633", "light sea green", "dodger blue", "#6CABE7", "gray")
    )

    #c <- read.table("cell_type_color.txt", sep="\t", check.names=F,  comment.char = "!", header=T)
    rownames(c) <- c$Cell.Type

    #cell type -> ident number
    cell_type_rank <- cbind(pbmc[[revised_head]][,1], Idents(pbmc))
    rownames(cell_type_rank) <- cell_type_rank[,1]
    cell_type_rank <- unique(cell_type_rank)
    colnames(cell_type_rank) <- c("Cell.Type", "Ident")

    print("Cell.Type")
    print(factor(x$Cell.Type))
    print("cell_type_rank")
    print(cell_type_rank)

     #ident -> color
    #level_order <- c("Cancer cells", "Adipocytes", "B cells", "Endothelial cells", "Normal epithelial cells", "Fibroblasts", "Macrophages", "Plasma cells", "T cells", "Pericytes")
    level_order <- c("Cancer cells", "Normal epithelial cells", "Stromal cells", "Immune cells")
    color_data <- merge(c, cell_type_rank)
    color_data <- color_data[order(as.numeric(as.character(color_data$Ident))),]
    color_flag <- as.vector(color_data$Color)
    color_label<- as.vector(color_data$Cell.Type)

    print("color_data")
    print(color_data)
    print("color_flag")
    print(color_flag)
    print("color_label")
    print(color_label)

    #color_flag <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")
    #color_label <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")

    #pdf("test.pdf", height=7, width=8)
    #ggplot(x, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(fill=factor(Cell.Type)), size=2, shape=21, stroke=0)
    p <- ggplot(x, aes(x= tSNE_1, y= tSNE_2)) + geom_point(aes(fill=factor(Cell.Type)), size=3, shape=21, stroke=0) + scale_fill_manual(values=color_flag, labels=color_label) + labs(fill="", colour="", stroke=0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.key = element_rect(fill = NA, colour = NA, size = 2)) + theme(text = element_text(size=18), plot.margin = margin(1, 1, 1, 1, "cm")) + guides(fill = guide_legend(override.aes = list(size=4)))
    filename1 <- paste(figdir, '/', fig_tsne, '.ggplot.broadclass.plot.pdf', sep="")
    pdf(filename1, width=8, height=7)
    p
    dev.off()

    font_size = 66
    p <- ggplot(x, aes(x= tSNE_1, y= tSNE_2)) + geom_point(aes(fill=factor(Cell.Type)), size=8, shape=21, stroke=0) + scale_fill_manual(values=color_flag, labels=color_label) + labs(fill="", colour="", stroke=0) + theme(plot.title=element_text(size=font_size,face="bold"), axis.text=element_text(size=font_size, face="bold"), axis.title=element_text(size=font_size,face="bold"), legend.text=element_text(size=font_size, face="bold"), legend.title=element_text(size=font_size, face="bold")) + theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.margin = margin(1, 1, 1, 1, "cm")) + guides(shape = guide_legend(override.aes = list(size = 15))) + theme(legend.key = element_blank())
    filename2 <- paste(figdir, '/', fig_tsne, '.ggplot.broadclass.png', sep="")
    png(filename2, width=2000, height=1200)
    p
    dev.off()
