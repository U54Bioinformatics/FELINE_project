rm(list=ls())
require(data.table);require(dplyr);require(ggplot2);require(tidyr);require(parallel);library(zinbwave);library( scRNAseq);library(biomaRt);library(umap);

args <- commandArgs(trailingOnly = TRUE)
cpu  <- 6
infile=args[1]
metafile=gsub('txt', 'metadata.txt', infile)
prefix=gsub('.txt', '', infile)
outfile_raw=paste0(prefix, '.zinbwave.raw.txt')
outfile_residuals=paste0(prefix, '.zinbwave.residuals.txt')
outfile_normalized=paste0(prefix, '.zinbwave.normalized.txt')
print("infile & metafile")
print(infile)
print(metafile)
pdf(paste0(prefix, '.zinbwave.umap.pdf'))
# load read data
print("load read data")
threads=cpu
setDTthreads(threads)
readcount_dd0 <- fread(infile, sep="\t", header=T)
readcount_dd0 <- data.table(readcount_dd0)
dim(readcount_dd0)
readcount_dd0[1:5,1:5]
#readcount_dd0 <- data.table(x)
#readcount_dd0 <- data.table(fread(infile, sep="\t", header=T))

# load meta and keep only these in meta file
#Cell.ID orig.ident      nCount_RNA      nFeature_RNA    Sample  Celltype_subtype        Celltype
#FEL011_M_AAAGTCCCAGCGAGTA       FEL011  5862    2691    FEL011P138_M    CAF-S1  Fibroblasts
#FEL011_M_AAAGTGAAGGGAGGGT       FEL011  4011    2218    FEL011P138_M    CAF-S1  Fibroblasts
#FEL011_M_AAATGGACACCTGATA       FEL011  3435    1846    FEL011P138_M    CAF-S1  Fibroblasts
meta_dd0 <- fread(metafile, sep="\t", header=T)
print("meta data")
dim(meta_dd0)
readcount_dd0 %>% dplyr::select(c('Gene.ID', meta_dd0$Cell.ID)) -> readcount_dd0
print("After filter with metadata")
dim(readcount_dd0)
readcount_dd0[1:5,1:5]

print("preprocess before normalization")
#remove unwanted genes
readcount_dd <- readcount_dd0 [! substr(Gene.ID,1,3)=="MT-"]

#Calc read depth per cell and ln read count
readcount_per_cell <- ( colSums(data.table( readcount_dd %>% dplyr::select(-"Gene.ID") )) )
lnreadcount_per_cell <- log( colSums(data.table( readcount_dd %>% dplyr::select(-"Gene.ID") )) )       #hist(readcount_per_cell, breaks=200)

# Initial gene list
genes<- as.factor(readcount_dd$Gene.ID)    #length(genes)

# select just the count data, not the gene name info, record cell ids and figure out the unique samples to add into the model as batch effects
rawcount_dd <- data.table(readcount_dd %>% dplyr::select(-Gene.ID))
rawcount_dd[1:5,1:5]
# get cell id names and put raw read count data into a matrix
cell_ids <- names(rawcount_dd)
rawcount_matrix <- as.matrix(rawcount_dd); rownames(rawcount_matrix)<-genes

#create single cell object and add cell specific info to SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = rawcount_matrix),
                            colData=data.frame(NREADS = readcount_per_cell, ln_NREADS = lnreadcount_per_cell,
                                               sample= as.factor(sapply(strsplit(cell_ids,"_"), `[`, 1))   ))

#filter rare genes :: Jeff may be able to optimize these settings:::  plot(sapply(1:40,function(x){sum(rowSums(assay(sce)>x)>10)} ))
filter <- rowSums(assay(sce)>1)>10 #table(filter)
table(filter)
sce <- sce[filter,]

#identify DE genes :: Jeff, do we want to limit the dataset to DE genes at this stage?
assay(sce) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)    #head(vars)
#filter differentially expressed genes if required
#sce_de <- sce[names(vars)[1:length(names(vars))],]#[1:500],]   # by default i keep all genes, but change to a specific n if wanting to filter for DE genes
#num_gene <- 2000
num_gene <- length(names(vars))
sce_de <- sce[names(vars)[1:num_gene],]
assayNames(sce_de)[1] <- "counts"
genes_sub <- rownames(sce_de)                    # length(names(vars));length(genes)

# get gene length and gc content :: this data could be collected for all genes and saved in csv, so that this step need not be run every time the analysis is required
#ensembl_gene_file='Ensembl_gene_info.txt'
#bm = ''
#if (!file.exists(ensembl_gene_file)) {
#    print("Get gene info from Ensembl with useMart")
    #You can also use asia.ensembl.org and uswest.ensembl.org for the host argument.
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
    bm <- getBM(attributes=c('hgnc_symbol', 'start_position','end_position', 'percentage_gene_gc_content'),
            filters = 'hgnc_symbol',
            values = rownames(sce_de),
            mart = mart)
    bm$length <- bm$end_position - bm$start_position
#    fwrite(bm, ensembl_gene_file, sep="\t", quote=F, col.names=T)
#} else {
#    print("Get gene info from local file: Ensembl_gene_info.txt")
#    bm = fread(ensembl_gene_file, sep="\t", header=T)
#}
len <- tapply(bm$length, bm$hgnc_symbol, mean)
len <- len[rownames(sce_de)]
gcc <- tapply(bm$percentage_gene_gc_content, bm$hgnc_symbol, mean)
gcc <- gcc[rownames(sce_de)]


#add gene specific info (and any further cell specific info) to SingleCellExperiment object
# database is not perfect and so genes with missing data are just assigned the expected values.
gene_covs <-   data.table(gccontent = as.vector(gcc), length = as.vector(len) )
gene_covs[!is.finite(length), length:=mean(len,na.rm=TRUE)]
gene_covs[!is.finite(gccontent), gccontent:=mean(gcc,na.rm=TRUE)]

rowData(sce_de) <- data.frame(gene_covs)  #colData(sce_de) <- data.frame(NREADS = readcount_per_cell, ln_NREADS = lnreadcount_per_cell)

# Covariate data
colData(sce_de)  # cell specific info
rowData(sce_de)  # gene specific info

# Test out some models :: Do this with small subset of genes and or cells for faster performance
#zz<-zinbFit(sce_de,   X="~ln_NREADS + sample",V="~gccontent + log(length)", K=2,epsilon=1000)
#zz2<-zinbFit(sce_de,   X="~ln_NREADS",V="~log(length)", K=2,epsilon=1000)
#zz3<-zinbFit(sce_de,   X="~NREADS",V="~log(length)", K=2,epsilon=1000)
#zz4<-zinbFit(sce_de,   X="~ln_NREADS",V="~(length)", K=2,epsilon=1000)
# dev_resid<- computeDevianceResiduals(zz, t(assays(sce_de)$counts), ignoreW = TRUE)
# model comparison to choose covariates
#zinbAIC(zz,t(assays(sce_de)$counts))
#zinbAIC(zz2,t(assays(sce_de)$counts))
#zinbAIC(zz3,t(assays(sce_de)$counts))
#zinbAIC(zz4,t(assays(sce_de)$counts))

#random set of cells
#set.seed(9)
#cell_random <- sample(seq_len(ncol(sce)), 200, replace=FALSE)
#print("random cells")
#sort(cell_random)
#sce_de_rn   <- sce_de[,cell_random]
#print("random cell obj")
#sce_de_rn

print("Zinbwave ormalization")
## Run the model to get normalised scores :: WARNING, this is getting slow with only moderately sized datasets.
start_time <- Sys.time()
start_time
#print("zinbFit & zinbwave")
#zinb_model <- zinbFit(sce_de_rn, K=2,  X="~ln_NREADS",V="~gccontent + log(length)", epsilon=1000)
#zinb_model <- zinbFit(sce_de_rn, K=2,  X="~ln_NREADS + sample",V="~gccontent + log(length)", epsilon=1000)
#print("use model to predict")
#zinb_bias_model <- zinbwave(sce_de, fitted_model = zinb_model, normalizedValues=TRUE, residuals = TRUE)
print("zinb fitted model")
#zinb_model
print("zinbwave multiple cpu")
zinb_bias_model <- zinbwave(sce_de, K=2, X="~ln_NREADS",V="~gccontent + log(length)", epsilon=1000, imputedValues=TRUE, normalizedValues=TRUE, residuals = TRUE, BPPARAM=MulticoreParam(cpu))
print("zinbwave single cpu")
#zinb_bias_model <- zinbwave(sce_de, K=2, X="~ln_NREADS + sample",V="~gccontent + log(length)", epsilon=1000, imputedValues=TRUE, normalizedValues=TRUE, residuals = TRUE)
print("zinbsurf")
#zinb_bias_model <- zinbsurf(sce_de, K=2, epsilon=1000, prop_fit = 0.1)
print("zinbwave return object")
zinb_bias_model
end_time <- Sys.time()
end_time
elaps<- end_time - start_time
elaps

# Outputs for the pipeline: People like Fred may want the raw residuals, the normalizedValues_residuals are used for further analysis and the raw counts may be wanted.
raw_residuals <- assays(zinb_bias_model)$residuals
normalizedValues_residuals <- assays(zinb_bias_model)$normalizedValues
raw_counts <- assays(zinb_bias_model)$counts

print("Write output")
out <- data.table(Gene.ID = names(vars)[1:num_gene])
fwrite(cbind(out, raw_counts), outfile_raw, sep="\t", quote=F, col.names=T, row.names=F)
fwrite(cbind(out, raw_residuals), outfile_residuals, sep="\t", quote=F, col.names=T, row.names=F)
fwrite(cbind(out, normalizedValues_residuals), outfile_normalized, sep="\t", quote=F, col.names=T, row.names=F)

print("umap")
##### Dimension reduction using raw and adjusted data
umap_data_raw <- umap( (t(assays(zinb_bias_model)$counts )))
umap_dd_raw <- data.table(umap_data_raw$layout);setnames(umap_dd_raw,old=paste0("V",1:umap_data_raw$ config$n_components),new=paste0("Dim",1:umap_data_raw$ config$n_components))
ggplot(umap_dd_raw,aes(Dim1, Dim2)) + geom_point(size=4) + # , shape=coverage
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()#+ facet_wrap(~bio)

umap_data <- umap( (t(assays(zinb_bias_model)$normalizedValues )))#umap_data <- umap( (t(assays(zinb_bias_model)$normalizedValues )) ,n_components = 2)
umap_dd <- data.table(umap_data$layout);setnames(umap_dd,old=paste0("V",1:umap_data$ config$n_components),new=paste0("Dim",1:umap_data$ config$n_components))
ggplot(umap_dd,aes(Dim1, Dim2)) + geom_point(size=4) + # , shape=coverage
  scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()#+ facet_wrap(~bio)


# what are these umap dimensions correlated with??
crr_dim1 <- cor(umap_data$layout[,1],t(assays(zinb_bias_model)$normalizedValues ))
crr_dim1_ord <- crr_dim1[order(-abs(crr_dim1))]
names(crr_dim1_ord) <- colnames(crr_dim1)[order(-abs(crr_dim1))]

crr_dim2 <- cor(umap_data$layout[,2],t(assays(zinb_bias_model)$normalizedValues ))
crr_dim2_ord <- crr_dim2[order(-abs(crr_dim2))]
names(crr_dim2_ord) <- colnames(crr_dim2)[order(-abs(crr_dim2))]

# first 30
(crr_dim1_ord)[1:30]
(crr_dim2_ord)[1:30]


# Gather final data and plot some overlays
ddFin<-data.frame(umap_dd,#Dim1=umap_data$layout[,1], Dim2=umap_data$layout[,2],
                  t(normalizedValues_residuals),
                  colData(zinb_bias_model))

ddFin_raw<-data.frame(umap_dd_raw,#Dim1=umap_data$layout[,1], Dim2=umap_data$layout[,2],
                      t(raw_counts),
                      colData(zinb_bias_model))

ggplot(ddFin_raw,aes(Dim1, Dim2,  col=ESR1)) + geom_point(size=4) + # , shape=coverage
  theme_classic()+ facet_wrap(~sample)

ggplot(ddFin,aes(Dim1, Dim2,  col=ESR1)) + geom_point(size=4) + # , shape=coverage
  theme_classic()+ facet_wrap(~sample)
dev.off()
print("zinbwave pipeline Done")
