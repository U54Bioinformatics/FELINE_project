args <- commandArgs(trailingOnly = TRUE)
prefix=args[1]
meta <- read.table(paste(prefix, "_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.txt", sep=""), header=T, sep="\t")
crosstable <- xtabs(~Sample+seurat_clusters, data=meta)
write.table(crosstable, file=paste(prefix, "_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.sample_matrix.txt", sep=""), sep="\t", quote=F)

