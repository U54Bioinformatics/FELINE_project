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

### Loading input data
# Founding clone must be cluster 1
# https://github.com/hdng/clonevol/issues/7 
args <- commandArgs(trailingOnly = TRUE)
prefix=args[1]
#sufix=args[2]
samples=args[2]
analysis=args[3]
samples_vector <- strsplit(samples, ",")[[1]]

#test
#prefix='test'
#sufix='short'

#infile=paste0(prefix, '.', sufix, 'txt')
infile=paste0(prefix, '.txt')
x <- read.table(infile, sep="\t", header=T)
## remove row if vaf is nan 
x <- x[!is.na(x[,ncol(x)]),]
print("x raw")
print(x[1:3, 1:ncol(x)])
vaf.col.names <- samples_vector 
sample.names  <- samples_vector
sample.groups <- samples_vector
names(sample.groups) <- sample.names

### testin data: AML1 clustering data
if(FALSE){
data(aml1)
x <- aml1$variants
print("x raw")
print(x[1:3, 1:ncol(x)])

### Output example data as short
write.table(x, "test.txt", sep="\t", quote=F, row.names=F)
x_short <- x[, c("cluster","gene", "is.driver", "P.ccf", "R.ccf", "P", "R")]
write.table(x_short, "test.short.txt", sep="\t", quote=F, row.names=F)

### Preparing clustering data
vaf.col.names <- grep('.vaf', colnames(x), value=T) 
sample.names <- gsub('.vaf', '', vaf.col.names)
x[, sample.names] <- x[, vaf.col.names]
vaf.col.names <- sample.names 
print("vaf.col.names, sample.names")
print(vaf.col.names)
print(sample.names)

### prepare sample grouping 
sample.groups <- c('P', 'R'); 
names(sample.groups) <- vaf.col.names 
print(vaf.col.names)
print(sample.names)
} #End of testing data

### setup the order of clusters to display in various plots (later) 
x <- x[order(x$cluster),]
print("x for ploting")
print(x[1:3, 1:ncol(x)])

### Choosing colors for the clones
# more colors
# fish = setCol(fish, col=c("#E86A10","#56A4AA", "#3A78C4","#F1AB00", "#2F992D", "#70738B", "grey"))
# fish = setCol(fish, col=c("#E86A10","#56A4AA", "#3A78C4","#F1AB00", "#2F992D", "grey"))
clone.colors <- c('#999793', '#8d4891', '#fe9536', '#d7352e', 'blue', "#F1AB00", "#2F992D", "#70738B")
#clone.colors <- NULL
#fishplot color 
clone.colors = c('gray', "#E86A10","#56A4AA", "#3A78C4","#F1AB00","#2F992D", '#8d4891', '#fe9536', '#d7352e', 'blue')

### Plots
#pdf(paste0(prefix, ".clonevol_figure.pdf"), width = 10, height = 7, useDingbats = FALSE, title='')
#layout(matrix(c(1,2,4,2,3,4), 2, 3, byrow = TRUE),
#   widths=c(1,1,2), heights=c(1,1))
if(TRUE){
pdf(paste0(prefix, ".clonevol_figure1.pdf"), width = 7, height = 7, useDingbats = FALSE, title='')
# Visualizing the variant clusters
plot.variant.clusters(x,
	cluster.col.name = 'cluster',
	show.cluster.size = FALSE,
	cluster.size.text.color = 'blue',
	vaf.col.names = vaf.col.names,
	vaf.limits = 100, sample.title.size = 20,
	violin = FALSE, box = FALSE,
	jitter = TRUE, jitter.shape = 1,
	jitter.color = clone.colors, 
	jitter.size = 3, jitter.alpha = 1, 
	jitter.center.method = 'median', 
	jitter.center.size = 1, 
	jitter.center.color = 'darkgray', 
	jitter.center.display.value = 'none', 
	highlight = 'is.driver', 
	highlight.shape = 21, 
	highlight.color = 'blue',
 	highlight.fill.color = 'green', 
	highlight.note.col.name = 'gene', 
	highlight.note.size = 2, 
	order.by.total.vaf = FALSE)
dev.off()
}

# Plotting pairwise VAFs or CCFs across samples
plot.pairwise(x, 
	col.names = vaf.col.names,
	out.prefix = paste0(prefix, '.clonevol_figure2'), 
	colors = clone.colors,
        xMinSmall=0, xMaxSmall=100,
        yMinSmall=0, yMaxSmall=100)

# Plotting mean/median of clusters across samples (cluster flow)
title_y = "Variant allele frequency (%)"
if ( analysis == 'ccf' ){
    title_y = "Cancer Cell Frequency (%)"
}
pdf(paste0(prefix, '.clonevol_figure3.pdf'), width = 7, height = 5, useDingbats = FALSE, title='')
plot.cluster.flow(x, 
	vaf.col.names = vaf.col.names, 
	sample.names  = sample.names,
        line.size = 2,
        x.title = '',
        y.title = title_y,
	colors = clone.colors)
dev.off()


### Inferring clonal evolution trees
y = ''
if (analysis == 'ccf'){
     print("Inferring clonal evolution trees using ccf ")
     y = infer.clonal.models(variants = x, 
	cluster.col.name = 'cluster', 
	ccf.col.names = vaf.col.names, 
	sample.groups = sample.groups, 
	cancer.initiation.model='monoclonal', 
	subclonal.test = 'bootstrap', 
	subclonal.test.model = 'non-parametric', 
	num.boots = 1000, 
	founding.cluster = 1, 
	cluster.center = 'mean', 
	ignore.clusters = NULL, 
	clone.colors = clone.colors, 
	min.cluster.vaf = 0.01, # min probability that CCF(clone) is non-negative 
	sum.p = 0.05, # alpha level in confidence interval estimate for CCF(clone) 
	alpha = 0.05)
} else {
     print("Inferring clonal evolution trees using vaf ")
     y = infer.clonal.models(variants = x,
        cluster.col.name = 'cluster', 
        vaf.col.names = vaf.col.names, 
        sample.groups = sample.groups, 
        cancer.initiation.model='monoclonal', #monoclonal, polyclonal 
        subclonal.test = 'bootstrap', 
        subclonal.test.model = 'non-parametric', 
        num.boots = 1000, 
        founding.cluster = 1, 
        cluster.center = 'mean', 
        ignore.clusters = NULL, 
        clone.colors = clone.colors,
        min.cluster.vaf = 0.01, # min probability that CCF(clone) is non-negative 
        sum.p = 0.05, # alpha level in confidence interval estimate for CCF(clone) 
        alpha = 0.05)
}

if(FALSE){
### Mapping driver events onto the trees
y <- transfer.events.to.consensus.trees(y, 
	x[x$is.driver,], 
	cluster.col.name = 'cluster', 
	event.col.name = 'gene')
}
y <- convert.consensus.tree.clone.to.branch(y, 
	branch.scale = 'sqrt')

### Fishplot
pdf(paste0(prefix, '.clonevol_figure4_fishplot.pdf'), width = 5, height = 3, useDingbats = FALSE)
par(mar=c(5, 4, 4, 2), xpd=TRUE)
fishplot_input <- generateFishplotInputs(y, rescale = TRUE)
fishplot_rdata = paste0(prefix, '.clonevol_figure4_fishplot.RData')
save(fishplot_input, sample.names, clone.colors, file=fishplot_rdata)

#clone.colors = c('gray', "#E86A10","#56A4AA", "#3A78C4","#F1AB00","#2F992D")
for (i in 1:length(fishplot_input$cell.fractions))
{
    
    fishplot_objects <- createFishObject(fishplot_input$cell.fractions[[i]], fishplot_input$parents[[i]])
    fishplot_layout = layoutClones(fishplot_objects)
    n_clone = length(fishplot_input$clonevol.clone.colors)
    print(n_clone) 
    fishplot_layout = setCol(fishplot_layout, clone.colors[1:n_clone])
    fishPlot(
        fishplot_layout,
        shape = "spline",
        cex.title = 0.7,
        vlines = seq(1, length(sample.names)),
        vlab = sample.names,
        pad.left = 0.2
    )
    legend("bottom", title="Cluster", c("1","2","3","4","5"),  fill=clone.colors, inset=c(0, -0.6), ncol=6, cex=1, bty="n")

}
dev.off()

### Plotting trees
pdf(paste0(prefix, '.clonevol_figure4_trees.pdf'), width = 3, height = 5, useDingbats = FALSE) 
plot.all.trees.clone.as.branch(y, 
	branch.width = 1, 
	node.size = 1.5, 
	node.label.size = 0.5) 
dev.off()

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
