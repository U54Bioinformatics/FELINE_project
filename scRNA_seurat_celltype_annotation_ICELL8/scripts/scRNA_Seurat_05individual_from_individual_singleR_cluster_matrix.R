args <- commandArgs(trailingOnly = TRUE)
prefix=args[1]
meta <- read.table(paste(prefix, "_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.txt", sep=""), header=T, sep="\t")
crosstable <- xtabs(~Encode_main_type+RNA_snn_res.0.8, data=meta)
write.table(crosstable, file=paste(prefix, "_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.SingleR_cluster_table_Blueprint.txt", sep=""), sep="\t", quote=F)

