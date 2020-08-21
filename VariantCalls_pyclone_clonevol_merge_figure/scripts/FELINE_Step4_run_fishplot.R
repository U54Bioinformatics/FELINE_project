#################################
#Match gene cnv with subclone tree:
#1. check cellular fraction of cnv of candidate genes in FACETS files: cf.em
#2. match with subclone clusters
#less -S ../VariantCalls_SNPs/FEL034_FACETSResults/facets.cval_200.ndepth_35.nbhd_250.nhet_15/FEL034_S.seg
#less -S ../VariantCalls_SNPs/FEL034_FACETSResults/facets.cval_200.ndepth_35.nbhd_250.nhet_15/FEL034_E.seg 
#In the output, cf, tcn, lcn are the initial estimates of cellular fraction,
#total and minor copy number estimates, and cf.em, tcn.em, lcn.em are the
#estimates by the mixture model optimized using the EM-algorithm. cf is used as
#initial values for the EM algorithm. For diploid normal segments (total copy=2,
#minor copy=1), we report cellular fraction as 1 (100\% normal). The logOR data
#for a segment are summarized using the square of expected log-odds-ratio
#(mafR column).
#https://github.com/mskcc/facets/issues/45
#################################
#no subclone model infered 
#https://github.com/hdng/clonevol/issues/11
#https://github.com/hdng/clonevol/issues/6
#################################

### Loading ClonEvol
library(clonevol)
library(fishplot)
library(tidyverse)

### Loading input data
# Founding clone must be cluster 1
# https://github.com/hdng/clonevol/issues/7 
args <- commandArgs(trailingOnly = TRUE)
prefix=args[1]
patient=args[2]

# load and select tumor size changes
clinic = read.table("FELINE_clinical_response.size_change.txt", sep="\t", header=T)
if (patient == "FEL036" || patient == "FEL046"){ 
   clinic %>% filter(orig.ident %in% patient, Day %in% c(31)) -> clinic
}else if (patient == "FEL044" || patient == "FEL045"){
   clinic %>% filter(orig.ident %in% patient, Day %in% c(91)) -> clinic
}else if (patient == "FEL031"){
   clinic %>% filter(orig.ident %in% patient, Day %in% c(151)) -> clinic
}else{
   clinic %>% filter(orig.ident %in% patient, Day %in% c(180)) -> clinic
}
print("tumor size change")
print(clinic)
print(clinic$Burden_change)
tumor_size_change = as.numeric(clinic$Burden_change)

## fishplot
### this will load fishplot_input, clone.colors and sample.names 
load(paste0(prefix, '.clonevol_figure4_fishplot.RData'))
pdf(paste0(prefix, '.clonevol_figure4_fishplot.pdf'), width = 4, height = 3, useDingbats = FALSE)
par(mar=c(1, 1, 2, 1), xpd=TRUE)

# sample.names
timepoints = c(0, 180)
if (patient == "FEL036" || patient == "FEL044" || patient == "FEL045" || patient == "FEL046"){
   sample.names = c("D0", "D14")
   timepoints   = c(0, 14)
}else{
   sample.names = c("D0", "D180")
   timepoints   = c(0, 180)
}

#clone.colors = c('gray', "#E86A10","#56A4AA", "#3A78C4","#F1AB00","#2F992D")
#for (i in 1:length(fishplot_input$cell.fractions))
for (i in 1:1)
{
    print(sample.names)
    print(fishplot_input$parents[[i]])
    print(timepoints)
    print(fishplot_input$cell.fractions[[i]])
    # adjust to tumor size
    data_m = fishplot_input$cell.fractions[[i]]
    data_m = data.frame(data_m)
    names(data_m) = c("Pre", "Post")
    if (tumor_size_change < 1){
       data_m$Post   = data_m$Post*tumor_size_change
    }else{
       data_m$Pre    = data_m$Pre/tumor_size_change
    }
    data_m        = as.matrix(data_m)
    print(data_m)
    #print(fishplot_input$cell.fractions[[i]]["FEL022_E"])
    #fishplot_objects <- createFishObject(fishplot_input$cell.fractions[[i]], fishplot_input$parents[[i]])
    fishplot_objects = createFishObject(data_m, fishplot_input$parents[[i]])
    fishplot_layout = layoutClones(fishplot_objects)
    n_clone = length(fishplot_input$clonevol.clone.colors)
    print(n_clone) 
    fishplot_layout = setCol(fishplot_layout, clone.colors[1:n_clone])
    fishPlot(
        fishplot_layout,
        shape = "spline",
        title = NULL,
        cex.title = 1.5,
        cex.vlab  = 1.5,
        vlines = seq(1, 2),
        vlab = sample.names,
        bg.type="gradient", 
        bg.col=c("bisque", "bisque", "bisque"),
        pad.left = 0.1
    )
    #legend("bottom", title="Cluster", c("1","2","3","4","5"),  fill=clone.colors, inset=c(0, -0.6), ncol=6, cex=1, bty="n")
}
dev.off()

### Plotting trees
if(FALSE){
pdf(paste0(prefix, '.clonevol_figure4_trees.pdf'), width = 3, height = 5, useDingbats = FALSE) 
plot.all.trees.clone.as.branch(y, 
	branch.width = 1, 
	node.size = 1.5, 
	node.label.size = 0.5) 
dev.off()
}

if(FALSE){ # plot one tree
y_tree = y$matched$merged.trees[[2]]
print(y_tree)
plot.tree.clone.as.branch(y_tree)
}
### Visualizing trees predicted by other tools (not working)
if(FALSE){
y = import.tree('tree.tsv', 'variants.tsv')
y = convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')
y <- transfer.events.to.consensus.trees(y, y$variants[y$variants$is.driver,], cluster.col.name = 'cluster', event.col.name = 'gene')
pdf('imported-tree.pdf', width=3, height=5, useDingbats=F) 
plot.all.trees.clone.as.branch(y, branch.width = 0.5, node.size = 1, node.label.size = 0.5) 
dev.off()
}
