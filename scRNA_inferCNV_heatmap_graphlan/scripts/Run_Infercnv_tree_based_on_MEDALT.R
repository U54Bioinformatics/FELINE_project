library(data.table)
library(fastcluster)
library(RColorBrewer)
library(ape)
library(ggtree)
library(tidyverse)
library(phytools)
library(treeio)

args <- commandArgs(trailingOnly = TRUE)
patient=args[1]

#cnv_file=paste0("./MEDALT_cnv/", patient, ".HMMi6.final.CNV.txt.CNV.txt")
cnv_file=paste0("./MEDALT_cnv/", patient, ".HMMi6.final.CNV.txt")
tree_file=paste0("./MEDALT_tree/", patient, ".tree")
tree_pdf=paste0("./MEDALT_tree/", patient, ".tree.pdf")
x= read.table(cnv_file)
x=t(x)
dend <- fastcluster::hclust(dist(x), method="ward.D2")
tre <- dend %>% as.phylo
# write tree file
write.tree(tre, file=tree_file)
# plot tree
pdf(tree_pdf)
ggtree(tre, continuous = TRUE, yscale = "trait") + scale_color_viridis_c() + theme_minimal()
dev.off()

