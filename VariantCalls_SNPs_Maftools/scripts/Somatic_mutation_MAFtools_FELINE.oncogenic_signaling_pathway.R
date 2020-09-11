
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
library(tidyverse)
clinic = read.table(laml.clin, header=T, sep="\t")
pdf(paste0(prefix, "_14_oncogenic_signaling_pathway.pdf"), width=10, height=7)
###########################Arm A#######################################
print("Oncogenic Signaling Pathways Arm A")
clinic %>% filter(ARM %in% c("A")) -> clinic_sub
clinic_sub_A = as.vector(clinic_sub$Tumor_Sample_Barcode)
lamlA = subsetMaf(maf = laml, tsb = clinic_sub_A, fields="Tumor_Sample_Barcode")
#lamlA = subsetMaf(maf = laml, query = "ARM == A")
lamlA
## Writes maf summary to an output file with basename laml.
write.mafSummary(maf = lamlA, basename = 'laml_ArmA')

# 14. Oncogenic Signaling Pathways
try({
p1 <- OncogenicPathways(maf = lamlA)
## Tumor suppressor genes are in red, and oncogenes are in blue font.
p2 <- PlotOncogenicPathways(maf = lamlA, pathways = "TGF-Beta")
p3 <- PlotOncogenicPathways(maf = lamlA, pathways = "TP53")
p4 <- PlotOncogenicPathways(maf = lamlA, pathways = "Cell_Cycle")
print(p1)
print(p2)
print(p3)
print(p4)
})
##################################################################

###########################Arm B#######################################
print("Oncogenic Signaling Pathways Arm B")
clinic %>% filter(ARM %in% c("B")) -> clinic_sub
clinic_sub_A = as.vector(clinic_sub$Tumor_Sample_Barcode)
lamlA = subsetMaf(maf = laml, tsb = clinic_sub_A, fields="Tumor_Sample_Barcode")
#lamlA = subsetMaf(maf = laml, query = "ARM == B")
lamlA
## Writes maf summary to an output file with basename laml.
write.mafSummary(maf = lamlA, basename = 'laml_ArmB')

# 14. Oncogenic Signaling Pathways
try({
p1 <- OncogenicPathways(maf = lamlA)
## Tumor suppressor genes are in red, and oncogenes are in blue font.
p2 <- PlotOncogenicPathways(maf = lamlA, pathways = "TGF-Beta")
p3 <- PlotOncogenicPathways(maf = lamlA, pathways = "TP53")
p4 <- PlotOncogenicPathways(maf = lamlA, pathways = "Cell_Cycle")
print(p1)
print(p2)
print(p3)
print(p4)
})
##################################################################

###########################Arm C#######################################
print("Oncogenic Signaling Pathways Arm C")
clinic %>% filter(ARM %in% c("C")) -> clinic_sub
clinic_sub_A = as.vector(clinic_sub$Tumor_Sample_Barcode)
lamlA = subsetMaf(maf = laml, tsb = clinic_sub_A, fields="Tumor_Sample_Barcode")
#lamlA = subsetMaf(maf = laml, query = "ARM == C")
lamlA
## Writes maf summary to an output file with basename laml.
write.mafSummary(maf = lamlA, basename = 'laml_ArmC')

# 14. Oncogenic Signaling Pathways
try({
p1 <- OncogenicPathways(maf = lamlA)
## Tumor suppressor genes are in red, and oncogenes are in blue font.
p2 <- PlotOncogenicPathways(maf = lamlA, pathways = "TGF-Beta")
p3 <- PlotOncogenicPathways(maf = lamlA, pathways = "TP53")
p4 <- PlotOncogenicPathways(maf = lamlA, pathways = "Cell_Cycle")
print(p1)
print(p2)
print(p3)
print(p4)
})
##################################################################
dev.off()


