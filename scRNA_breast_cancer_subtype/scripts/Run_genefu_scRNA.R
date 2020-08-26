library(genefu)
library(xtable)
library(rmeta)
library(Biobase)
library(caret)
library(data.table)
library(matrixStats)

#Rscript Run_genefu_scRNA.R CPM.txt Celltype.txt
#
#Input file:
#CPM.txt:
#Gene.ID S_AAACCCAAGTATGACA       S_AAACCCACATTGAAGA
#LINC00115           0.00             0.00
#LINC02593           0.00             0.00

#Celltype.txt
#Cell.ID Sample  Celltype1       Celltype2
#M_AAACGAACACAAGTGG       M    Cancer cells    Epithelial cells
#M_AAACGCTGTTAAGAAC       M    Cancer cells    Epithelial cells
#M_AAAGGTATCAGCCTCT       M    Macrophages     Macrophages
#M_AAAGTCCCAGCGAGTA       M    CAF-S1  Fibroblasts


args <- commandArgs(trailingOnly = TRUE)
cpm_file=args[1]
celltype=args[2]
gene_id=args[3]
#cpm_file="FEL011040_10x_FEL011_gene_symbols.CPM.txt"
#celltype="FEL011046.cell_type.txt"
#gene_id="FEL011046.gene_id.txt"

#prefix="FEL011040_10x_FEL024_gene_symbols"
prefix=gsub('.CPM.txt', '', cpm_file)
prefix
sum_file <- paste0(prefix, ".genefu.sum.txt")
sum_ct1_file <- paste0(prefix, ".genefu.sum.celltype1.txt")
sum_ct2_file <- paste0(prefix, ".genefu.sum.celltype2.txt")
sum_sample_file <- paste0(prefix, ".genefu.sample_sum.txt")
pc_sample_file  <- paste0(prefix, ".genefu.sample_pc.txt")
sum_file1 <- paste0(prefix, ".genefu_filter.sum.txt")
sum_ct1_file1 <- paste0(prefix, ".genefu_filter.sum.celltype1.txt")
sum_ct2_file1 <- paste0(prefix, ".genefu_filter.sum.celltype2.txt")
sum_sample_file1 <- paste0(prefix, ".genefu_filter.sample_sum.txt")
pc_sample_file1  <- paste0(prefix, ".genefu_filter.sample_pc.txt")
cpm <- fread(cpm_file, header=T, sep="\t")
cell_type <- fread(celltype, sep="\t", header=T)


#load gene id
gene_id_raw = fread(gene_id, sep="\t")
names(gene_id_raw) = c("Gene.ID_1", "Gene.ID_2", "Gene.ID")

#save gene id
cpm_with_id = merge(cpm, gene_id_raw, on="Gene.ID")
genes <- data.table(Gene.Symbol=cpm_with_id$Gene.ID, EntrezGene.ID=cpm_with_id$Gene.ID_2)
genes <- data.frame(genes)
rownames(genes) = genes$Gene.Symbol
cpm <- cpm_with_id
cpm$Gene.ID <- NULL
cpm$Gene.ID_1 <- NULL
cpm$Gene.ID_2 <- NULL

#prepare log1(cpm) matrix
cpm_m <- as.matrix(cpm)
rownames(cpm_m) <- genes$Gene.Symbol
cpm_m_log1 <- t(log1p(cpm_m))

#run prediction
#SubtypePredictions<-molecular.subtyping(sbt.model = "pam50", data = cpm_m_log1, anno=genes, do.mapping = F)
#SubtypePredictions<-molecular.subtyping(sbt.model = "scmod1", data = cpm_m_log1, anno=genes, do.mapping = T)
SubtypePredictions<-molecular.subtyping(sbt.model = "scmod2", data = cpm_m_log1, anno=genes, do.mapping = T)
max_p <- data.table(Max_p=rowMaxs(SubtypePredictions$subtype.proba))

#output table
subtype <- data.table(Cell.ID=names(SubtypePredictions$subtype), Subtype=as.vector(SubtypePredictions$subtype))
x <- cbind(subtype, SubtypePredictions$subtype.proba, max_p, SubtypePredictions$subtype.crisp)
#pam50
#names(x) <- c("Cell.ID", "Subtype", "Basal_p", "Her2_p", "LumA_p", "LumB_p", "Normal_p", "Max_p", "Basal", "Her2", "LumA", "LumB", "Normal")
#x_filter <- x[x$Max_p >= 0.7]
#scmod1
names(x) <- c("Cell.ID", "Subtype", "ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif", "Max_p", "Basal", "Her2", "LumB", "LumA")
x_filter <- x
#############################raw####################################################
print("dim x")
print(dim(x))
print(dim(cell_type))
x_cell_type <- merge(x, cell_type, by="Cell.ID")
fwrite(x_cell_type, sum_file, sep="\t", quote=F, row.name=F)
print(dim(x_cell_type))

#output sum, cell type1
x_ct1_sum <- as.data.frame(table(x_cell_type$Subtype, x_cell_type$Celltype1))
names(x_ct1_sum) <- c("Subtype", "Celltype1", "Freq")
x_ct1_sum <- dcast(x_ct1_sum, Celltype1~Subtype)
fwrite(x_ct1_sum, sum_ct1_file, sep="\t", quote=F, row.name=F)

#output sum, cell type2
x_ct2_sum <- as.data.frame(table(x_cell_type$Subtype, x_cell_type$Celltype2))
names(x_ct2_sum) <- c("Subtype", "Celltype2", "Freq")
x_ct2_sum <- dcast(x_ct2_sum, Celltype2~Subtype)
fwrite(x_ct2_sum, sum_ct2_file, sep="\t", quote=F, row.name=F)

#output sample sum
x_cell_type_tumor <- x_cell_type[Celltype1=="Cancer cells"]
x_sample_sum <- as.data.frame(table(x_cell_type_tumor$Sample, x_cell_type_tumor$Subtype))
names(x_sample_sum) <- c("Sample", "Subtype", "Freq")
x_sample_sum <- dcast(x_sample_sum, Sample~Subtype)
x_sample_total <- data.table(Total=rowSums(as.matrix(x_sample_sum[,-1])))
x_sample_sum_t <- cbind(x_sample_sum, x_sample_total)
fwrite(x_sample_sum_t, sum_sample_file, sep="\t", quote=F, row.name=F)
#percent
x_sample_pc <- t(apply(x_sample_sum[,2:ncol(x_sample_sum)], 1, function(x){100*x/sum(x)}))
x_sample <- data.table(Sample=x_sample_sum_t$Sample)
x_sample_sum_pc <- cbind(x_sample, x_sample_pc, x_sample_total)
fwrite(x_sample_sum_pc, pc_sample_file, sep="\t", quote=F, row.name=F)
#############################filter with 0.7###################################################
if(FALSE){
x_cell_type <- merge(x_filter, cell_type, by="Cell.ID")
fwrite(x_cell_type, sum_file1, sep="\t", quote=F, row.name=F)

#output sum, cell type1
x_ct1_sum <- as.data.frame(table(x_cell_type$Subtype, x_cell_type$Celltype1))
names(x_ct1_sum) <- c("Subtype", "Celltype1", "Freq")
x_ct1_sum <- dcast(x_ct1_sum, Celltype1~Subtype)
fwrite(x_ct1_sum, sum_ct1_file1, sep="\t", quote=F, row.name=F)

#output sum, cell type2
x_ct2_sum <- as.data.frame(table(x_cell_type$Subtype, x_cell_type$Celltype2))
names(x_ct2_sum) <- c("Subtype", "Celltype2", "Freq")
x_ct2_sum <- dcast(x_ct2_sum, Celltype2~Subtype)
fwrite(x_ct2_sum, sum_ct2_file1, sep="\t", quote=F, row.name=F)

#output sample sum
x_cell_type_tumor <- x_cell_type[Celltype1=="Cancer cells"]
x_sample_sum <- as.data.frame(table(x_cell_type_tumor$Sample, x_cell_type_tumor$Subtype))
names(x_sample_sum) <- c("Sample", "Subtype", "Freq")
x_sample_sum <- dcast(x_sample_sum, Sample~Subtype)
x_sample_total <- data.table(Total=rowSums(as.matrix(x_sample_sum[,-1])))
x_sample_sum_t <- cbind(x_sample_sum, x_sample_total)
fwrite(x_sample_sum_t, sum_sample_file1, sep="\t", quote=F, row.name=F)
#percent
x_sample_pc <- t(apply(x_sample_sum[,2:ncol(x_sample_sum)], 1, function(x){100*x/sum(x)}))
x_sample <- data.table(Sample=x_sample_sum_t$Sample)
x_sample_sum_pc <- cbind(x_sample, x_sample_pc, x_sample_total)
fwrite(x_sample_sum_pc, pc_sample_file1, sep="\t", quote=F, row.name=F)
}
#####################################################################################
