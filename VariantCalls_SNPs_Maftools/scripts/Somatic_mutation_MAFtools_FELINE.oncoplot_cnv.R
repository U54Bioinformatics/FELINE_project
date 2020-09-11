
if (FALSE){
  # update bioconductor
  #if (!requireNamespace("BiocManager", quietly = TRUE))
  #    install.packages("BiocManager")
  #BiocManager::install(version = "3.10")
  # install maftools
  if (!require("BiocManager"))
    install.packages("BiocManager")
  BiocManager::install("maftools")
  # insall BSgenome.Hsapiens.UCSC.hg19
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19") 
}

library(maftools)
library(dplyr)
# Requires BSgenome object
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
args <- commandArgs(trailingOnly = TRUE)
prefix=args[1]
laml.maf=args[2]
laml.clin=args[3]
laml.cnv=args[4]
#prefix="FELINE"

# Load input files
## path to TCGA LAML MAF file
#laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
## clinical information containing survival information and histology. This is optional
#laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
# Load FELINE files
#laml.maf = "FELINE.vep.maf"
#laml.clin = "FELINE.clinical.txt"
#laml.cnv  = "FELINE_FACETS_copy_number_gene.cnv4maftools.txt"
#laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
laml = read.maf(maf = laml.maf, clinicalData = laml.clin, cnTable = laml.cnv)
laml
## Shows sample summry.
getSampleSummary(laml)
## Shows gene summary.
#getGeneSummary(laml)
#shows clinical data associated with samples
#getClinicalData(laml)
## Shows all fields in MAF
#getFields(laml)
###################filter data
#these have to removed if analyzing copy numbers because these samples have poor tumor purity to estimate cnv
#blacklist <- c("FEL019_S", "FEL019_E", "FEL029_S", "FEL029_E", "FEL037_S", "FEL037_E", "FEL014_E", "FEL021_E", "FEL023_E", "FEL027_S", "FEL046_S")
#13,36,43,25,31
#blacklist <- c("FEL019_S", "FEL019_E", "FEL029_S", "FEL029_E", "FEL037_S", "FEL037_E", "FEL043_S", "FEL043_E", "FEL014_E", "FEL021_E", "FEL023_E", "FEL027_S", "FEL031_E", "FEL046_S")
#laml.sample <- getSampleSummary(laml)
#samplelist <- laml.sample$Tumor_Sample_Barcode 
#samplelist <- samplelist[!samplelist %in% blacklist]

################################################################
clinic = read.table(laml.clin, header=T, sep="\t")
clinic %>% arrange(ARM, Response, Tumor_Sample_Barcode) -> clinic
sample_level <- as.vector(clinic$Tumor_Sample_Barcode)
clinic
################################################################
try({
pdf(paste0(prefix, "_02_oncoplots.pdf"), width=10, height=10)
par(mar=c(6,6,6,6))
oncoplot(maf = laml, clinicalFeatures = c("ARM", 'Response'), altered=TRUE, top = 25, fontSize=0.75, barcode_mar  = 2, gene_mar = 6, showTumorSampleBarcodes = TRUE)
oncoplot(maf = laml, clinicalFeatures = c("ARM", 'Response'), altered=FALSE, top = 25, fontSize=0.75, barcode_mar  = 2, gene_mar = 6)
# Enrich in NR
oncostrip(maf = laml, genes = c('PIK3CA', 'TP53', 'CDK6', 'CCND1', 'CCNE2', 'AKT1', 'AKT3', 'MYC', 'AURKA', 'ERBB2', 'ERBB4', 'FGFR1', 'FGFR2', 'FGF1', 'FGF4', 'FGF12', 'ESR1', 'NF1', 'RB1', 'FAT1'),
 showTumorSampleBarcodes=T, keepGeneOrder=T, clinicalFeatures = c("ARM", 'Response'), barcode_mar = 6, removeNonMutated=F, legendFontSize=1.5, annotationFontSize=1.5, sampleOrder = sample_level)
# oncogeic pathway
oncostrip(maf = laml, genes = c('RB1', 'CCND1', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CDK6', 'CCND3', 'CDK4', 'CDKN1A', 'CDKN1B', 'E2F3', 'CCND2', 'CCNE1', 'CDK2', 'TP53', 'ATM', 'MDM4', 'CHEK2', 'MDM2', 'SMAD4', 'SMAD2', 'SMAD3', 'ACVR2A', 'TGFBR2', 'ACVR1B'),
 showTumorSampleBarcodes=T, keepGeneOrder=T, clinicalFeatures = c("ARM", 'Response'), barcode_mar = 6, removeNonMutated=F, legendFontSize=1.5, annotationFontSize=1.5, sampleOrder = sample_level)

})
##################################################################
dev.off()


