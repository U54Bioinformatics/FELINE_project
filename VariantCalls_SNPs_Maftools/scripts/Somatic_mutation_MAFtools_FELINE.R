
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
getGeneSummary(laml)
#shows clinical data associated with samples
getClinicalData(laml)
## Shows all fields in MAF
getFields(laml)
###################filter data
#these have to removed if analyzing copy numbers because these samples have poor tumor purity to estimate cnv
#blacklist <- c("FEL019_S", "FEL019_E", "FEL029_S", "FEL029_E", "FEL037_S", "FEL037_E", "FEL014_E", "FEL021_E", "FEL023_E", "FEL027_S", "FEL046_S")
#13,36,43,25,31
#blacklist <- c("FEL019_S", "FEL019_E", "FEL029_S", "FEL029_E", "FEL037_S", "FEL037_E", "FEL043_S", "FEL043_E", "FEL014_E", "FEL021_E", "FEL023_E", "FEL027_S", "FEL031_E", "FEL046_S")
laml.sample <- getSampleSummary(laml)
samplelist <- laml.sample$Tumor_Sample_Barcode 
#samplelist <- samplelist[!samplelist %in% blacklist]
laml = subsetMaf(maf = laml, tsb = samplelist, fields="Tumor_Sample_Barcode")
###################
laml
## Shows sample summry.
getSampleSummary(laml)
## Shows gene summary.
getGeneSummary(laml)
#shows clinical data associated with samples
getClinicalData(laml)
## Shows all fields in MAF
getFields(laml)
## Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')

# 1. Plotting MAF summary.
try({
pdf(paste0(prefix, "_01_maf_summary.pdf"), width=12, height=10)
p1 <- plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', top = 25, dashboard = TRUE, titvRaw = FALSE, textSize = 1)
p2 <- plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = FALSE, titvRaw = FALSE)
print(p1)
print(p2)
dev.off()
})

# 2. Drawing oncoplots
sample_level = c("FEL014_S", "FEL014_E", "FEL029_S", "FEL029_E", "FEL031_S", "FEL031_E", "FEL037_S", "FEL037_E", "FEL021_S", "FEL021_E", "FEL022_S", "FEL022_E", "FEL023_S", "FEL023_E", "FEL025_S", "FEL025_E", "FEL027_S", "FEL027_E", "FEL028_S", "FEL028_E", "FEL045_S", "FEL045_M", "FEL046_S", "FEL046_M", "FEL013_S", "FEL013_E", "FEL015_S", "FEL015_E", "FEL019_S", "FEL019_E", "FEL020_S", "FEL020_E", "FEL026_S", "FEL026_E", "FEL034_S", "FEL034_E", "FEL035_S", "FEL035_E", "FEL036_S", "FEL036_M", "FEL039_S", "FEL039_E", "FEL041_S", "FEL041_E", "FEL043_S", "FEL043_E", "FEL044_S", "FEL044_M")
try({
pdf(paste0(prefix, "_02_oncoplots.pdf"), width=10, height=10)
par(mar=c(6,6,6,6))
oncoplot(maf = laml, clinicalFeatures = c("ARM", 'Response'), altered=TRUE, top = 25, fontSize=0.75, barcode_mar  = 2, gene_mar = 6, showTumorSampleBarcodes = TRUE)
oncoplot(maf = laml, clinicalFeatures = c("ARM", 'Response'), altered=FALSE, top = 25, fontSize=0.75, barcode_mar  = 2, gene_mar = 6)
## visualize any set of genes using oncostrip function
# resistant
oncostrip(maf = laml, genes = c('PIK3CA','TP53', 'ESR1', 'ERBB2', 'ERBB4', 'NF1', 'RB1', 'FAT1', 'AURKA', 'MDM2', 'CCND1', 'CDK2', 'CDK4', 'CDK6', 'CCNE1', 'CCNE2', 'AKT1', 'AKT2', 'AKT3', 'MYC', 'FGFR1', 'FGFR2', 'FGF1', 'FGF2', 'FGF3', 'FGF4'),
 showTumorSampleBarcodes=T, keepGeneOrder=T, clinicalFeatures = c("ARM", 'Response', 'Ki67_Response'), barcode_mar = 6, sampleOrder = sample_level, removeNonMutated=F, legendFontSize=1.5, annotationFontSize=1.5)
# GF
oncostrip(maf = laml, genes = c('ATF3', 'BTG2', 'CALR', 'FGFR1', 'FGF12', 'FGF14', 'IGF2R', 'IGFN1', 'ILF3', 'JUNB', 'KLF10', 'SEPHS1', 'TNFRSF12A', 'TNFSF13B'),
 showTumorSampleBarcodes=T, keepGeneOrder=T, clinicalFeatures = c("ARM", 'Response', 'Ki67_Response'), barcode_mar = 6, sampleOrder = sample_level, removeNonMutated=F, legendFontSize=1.5, annotationFontSize=1.5)
# Enrich in NR
oncostrip(maf = laml, genes = c('PIK3CA', 'TP53', 'CDK6', 'CCND1', 'CCNE2', 'AKT1', 'AKT3', 'MYC', 'AURKA', 'ERBB2', 'ERBB4', 'FGFR1', 'FGFR2', 'FGF1', 'FGF4', 'FGF12', 'ESR1', 'NF1', 'RB1', 'FAT1'),
 showTumorSampleBarcodes=T, keepGeneOrder=T, clinicalFeatures = c("ARM", 'Response', 'Ki67_Response'), barcode_mar = 6, sampleOrder = sample_level, removeNonMutated=F, legendFontSize=1.5, annotationFontSize=1.5)
# FGF
oncostrip(maf = laml, genes = c('FGFR1','FGFR2', 'FGFR3', 'FGFR4', 'FGF1', 'FGF2', 'FGF3', 'FGF4', 'FGF5', 'FGF6', 'FGF7', 'FGF8', 'FGF9', 'FGF10', 'FGF11', 'FGF12', 'FGF13', 'FGF14', 'FGF16', 'FGF17', 'FGF18', 'FGF19', 'FGF20', 'FGF21', 'FGF22', 'FGF23'),
 showTumorSampleBarcodes=T, keepGeneOrder=T, clinicalFeatures = c("ARM", 'Response', 'Ki67_Response'), barcode_mar = 6, sampleOrder = sample_level, removeNonMutated=F, legendFontSize=1.5, annotationFontSize=1.5)
# Top mutated genes
oncostrip(maf = laml, genes = c('PIK3CA', 'TP53', 'TTN', 'CDC27', 'MAP3K1', 'MUC16', 'PER3', 'AKAP9', 'CACNA1B', 'CCDC108', 'DDX18', 'LRP1B', 'MUC12', 'NCOR1', 'PDS5B', 'RADIL', "RBMXL3", "TLN2", "HGFAC", "THEMIS2"),
 showTumorSampleBarcodes=T, keepGeneOrder=T, clinicalFeatures = c("ARM", 'Response', 'Ki67_Response'), barcode_mar = 6, sampleOrder = sample_level, removeNonMutated=F, legendFontSize=1.5, annotationFontSize=1.5)
# Done
dev.off()

})

# 3. Transition and Transversions.
try({
pdf(paste0(prefix, "_03_titv.pdf"), width=10, height=7)
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
## plot titv summary
p1 <- plotTiTv(res = laml.titv)
## ti/tv ratio
print(summary(laml.titv$TiTv.fractions$Ti/laml.titv$TiTv.fractions$Tv))
p2 <- plot(density(laml.titv$TiTv.fractions$Ti/laml.titv$TiTv.fractions$Tv), main="Ti/Tv ratio", xlab="Ti/Tv ratio")
print(p1)
print(p2)
dev.off()
})

# 4. Lollipop plots for amino acid changes and Labelling points
try({
pdf(paste0(prefix, "_04_lollipop_genes.pdf"), width=10, height=7)
p1 <- lollipopPlot(maf = laml, gene = 'PIK3CA', AACol = 'HGVSp_Short', showMutationRate = TRUE)
p2 <- lollipopPlot(maf = laml, gene = 'PER3', AACol = 'HGVSp_Short', showMutationRate = TRUE)
p3 <- lollipopPlot(maf = laml, gene = 'TTN', AACol = 'HGVSp_Short', showMutationRate = TRUE)
print(p1)
print(p2)
print(p3)
dev.off()
})

# 5. rainfallPlot
try({
pdf(paste0(prefix, "_05_rainfall.pdf"), width=10, height=7)
#p1 <- rainfallPlot(maf = laml, detectChangePoints = TRUE, pointSize = 0.6)
#print(p1)
dev.off()
})

# 6. Compare mutation load against TCGA cohorts
try({
pdf(paste0(prefix, "_06_mutation_load.pdf"), width=10, height=7)
try(tcgaCompare(maf = laml, cohortName = 'FELINE', logscale=TRUE))
dev.off()
})

# 7. gene vaf and Genecloud
try({
pdf(paste0(prefix, "_07_gene_vaf_clound.pdf"), width=10, height=7)
try(plotVaf(maf = laml, vafCol = 'i_TumorVAF_WU'))
try(geneCloud(input = laml, minMut = 3))
dev.off()
})

# 8. Somatic Interactions
##exclusive/co-occurance event analysis on top 10 mutated genes. 
try({
pdf(paste0(prefix, "_08_somatic_interaction.pdf"), width=10, height=7)
par(mar=c(6,6,4,2))
p1 <- somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.65)
print(p1)
dev.off()
})

# 9. Detecting cancer driver genes based on positional clustering
try({
pdf(paste0(prefix, "_09_clustered_driver_genes.pdf"), width=10, height=7)
laml.sig = oncodrive(maf = laml, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
head(laml.sig)
p1 <- plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)
print(p1)
dev.off()
})

# 10. Adding and summarizing pfam domains
try({
pdf(paste0(prefix, "_10_pfam.pdf"), width=10, height=7)
laml.pfam = pfamDomains(maf = laml, AACol = 'HGVSp_Short', top = 10)
dev.off()
## Protein summary (Printing first 7 columns for display convenience)
laml.pfam$proteinSummary[,1:7, with = FALSE]
## Domain summary (Printing first 3 columns for display convenience)
laml.pfam$domainSummary[,1:3, with = FALSE]
})

# 11. Pan-Cancer comparison
# Lawrence et al performed MutSigCV analysis on 21 cancer cohorts and 
# identified over 200 genes to be significantly mutated which consists of
# previously un-subscribed novel genes 9. Their results show only few genes
# are mutated in multiple cohort while many of them are tissue/cohort specific. 
# We can compare mutSig results against this pan-can list of significantly 
# mutated genes to see genes specifically mutated in given cohort. 
# This function requires MutSigCV results (usually named sig_genes.txt) as an input.
# This function takes MutSig results and compares them against panCancer cohort (~5000 tumor samples from 21 cancer types). 
# This analysis can reveal novel genes exclusively mutated in input cohort.
# MutsigCV results for TCGA-AML
#try({
#laml.mutsig <- system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
#pdf(paste0(prefix, "_11_pan_cancer_comparison.pdf"), width=10, height=7)
#p1 <- pancanComparison(mutsigResults = laml.mutsig, qval = 0.1, cohortName = 'FELINE', inputSampleSize = 200, label = 1)
#print(p1)
#dev.off()
#})

# 12. Clinical enrichment analysis
try({
fab.ce = clinicalEnrichment(maf = laml, clinicalFeature = 'Response')
## Results are returned as a list. Significant associations p-value < 0.05
fab.ce$groupwise_comparision[p_value < 0.05]
write.table(fab.ce$groupwise_comparision[p_value < 0.05], paste0(prefix, "_12_clinical_enrich_response.txt"), quote=FALSE, sep="\t")
pdf(paste0(prefix, "_12_clinical_enrich_response.pdf"), width=10, height=7)
p1 <- plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05)
print(p1)
dev.off()
})

try({
fab.ce = clinicalEnrichment(maf = laml, clinicalFeature = 'ARM')
## Results are returned as a list. Significant associations p-value < 0.05
fab.ce$groupwise_comparision[p_value < 0.05]
write.table(fab.ce$groupwise_comparision[p_value < 0.05], paste0(prefix, "_12_clinical_enrich_arm.txt"), quote=FALSE, sep="\t")
pdf(paste0(prefix, "_12_clinical_enrich_arm.pdf"), width=10, height=7)
p1 <- plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05)
print(p1)
dev.off()
})


# 13. Drug-Gene Interactions
try({
dgi = drugInteractions(maf = laml, fontSize = 0.75)
pik3ca.dgi = drugInteractions(genes = "PIK3CA", drugs = TRUE)
per3.dgi = drugInteractions(genes = "PER3", drugs = TRUE)
ttn.dgi = drugInteractions(genes = "TTN", drugs = TRUE)
## Printing selected columns.
pik3ca.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]
per3.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]
ttn.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]
})

# 14. Oncogenic Signaling Pathways
try({
pdf(paste0(prefix, "_14_oncogenic_signaling_pathway.pdf"), width=10, height=7)
p1 <- OncogenicPathways(maf = laml)
## Tumor suppressor genes are in red, and oncogenes are in blue font.
p2 <- PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS")
p3 <- PlotOncogenicPathways(maf = laml, pathways = "Hippo")
print(p1)
print(p2)
print(p3)
dev.off()
})

# 15. APOBEC Enrichment estimation
try({
laml.tnm = trinucleotideMatrix(maf = laml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
pdf(paste0(prefix, "_15_apobec_enrich.pdf"), width=10, height=7)
p1 <- plotApobecDiff(tnm = laml.tnm, maf = laml, pVal = 0.2)
print(p1)
dev.off()
})


