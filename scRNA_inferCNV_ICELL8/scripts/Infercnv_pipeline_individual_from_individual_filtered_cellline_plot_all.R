library(infercnv)
args = commandArgs(trailingOnly=TRUE)
#Sys.setenv("DISPLAY"=":0")
#options(bitmapType='cairo')
# create the infercnv object
sample=args[1]
platform=args[2]
down=args[3]
control=args[4]
cell=as.numeric(args[5])
prefix=paste(sample, platform, sep="_")
out=paste(prefix, "_counts.", down, ".output_dir_plot_", cell, "cells", sep="")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste(prefix, "counts.matrix", sep="_"),
                                     annotations_file=paste(prefix, "_counts.anno.", down, ".txt", sep=""),
                                     delim="\t",
                                     gene_order_file=paste(prefix, "counts.gene_pos.txt", sep="_"),
                                     ref_group_names=c("Immune"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                              cutoff=0,
                              min_cells_per_gene = cell,
                              out_dir=out, 
                              cluster_by_groups=T, 
                              denoise=T,
                              HMM=T,
                              debug = FALSE,
                              num_threads = 16,
                              plot_steps = F,
                              no_plot = T,
                              plot_probabilities = F,
                              png_res = 300,
                              analysis_mode="subclusters"
 )

#output obj 1
obj_1="14_HMM_pred.repr_intensitiesHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj"
out_1="infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities"
in_obj=paste0(out, "/", obj_1)
out_dir=out
out_file=out_1
infercnv_obj = readRDS(in_obj)
plot_cnv(infercnv_obj,
         out_dir=out_dir,
         output_filename=out_file,
         x.range="auto",
         x.center=1,
         title="",
         output_format=NA,
         color_safe_pal = FALSE)

#output obj 2
obj_2="run.final.infercnv_obj"
out_2="infercnv"
in_obj=paste0(out, "/", obj_2)
out_dir=out
out_file=out_2
infercnv_obj = readRDS(in_obj)
plot_cnv(infercnv_obj,
         out_dir=out_dir,
         output_filename=out_file,
         x.range="auto",
         x.center=1,
         title="",
         output_format=NA,
         color_safe_pal = FALSE)

#output obj 3
obj_3="preliminary.infercnv_obj"
out_3="infercnv.preliminary"
in_obj=paste0(out, "/", obj_3)
out_dir=out
out_file=out_3
infercnv_obj = readRDS(in_obj)
plot_cnv(infercnv_obj,
         out_dir=out_dir,
         output_filename=out_file,
         x.range="auto",
         x.center=1,
         title="",
         output_format=NA,
         color_safe_pal = FALSE)

