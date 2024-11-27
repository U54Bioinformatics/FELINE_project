# Breast cancer cells communicate with macrophages to prevent T cell activation during development of cell cycle therapy resistance 

This repository contains code accompanying the manuscript: "Breast cancer cells communicate with macrophages to prevent T cell activation during development of cell cycle therapy resistance". It provides code used to analyze the breast tumor composition, explore cellular phenotypes linked to treatment resistance and decipher tumor-wide communications from phenotypically diverse populations of cancer and non-cancer cells within the tumor microenvironment (TME). Serial ER+ breast tumor biopsies from of post-menopausal women on the FELINE clinical trial (clinicaltrials.gov # NCT02712723) were profiled using single cell RNA-seq (scRNAseq). Tumor samples from two cohorts of patients were independently examined, with the first discovery cohort used to identify tumor features promoting treatment resistance. The second cohort was used to validate discoveries and critical predictions. 


Related publications (Please consult and cite accordingly when using this repository):

Griffiths, J. I., Chen, J., Cosgrove, P. A., O’Dea, A., Sharma, P., Ma, C., ... & Bild, A. H. (2021). Serial single-cell genomics reveals convergent subclonal evolution of resistance as patients with early-stage breast cancer progress on endocrine plus CDK4/6 therapy. Nature cancer, 2(6), 658-671.https://www.nature.com/articles/s43018-021-00215-7

Griffiths, J. I., Wallet, P., Pflieger, L. T., Stenehjem, D., Liu, X., Cosgrove, P. A., ... & Bild, A. H. (2020). Circulating immune cell phenotype dynamics reflect the strength of tumor–immune cell interactions in patients during immunotherapy. Proceedings of the National Academy of Sciences, 117(27), 16072-16082. https://www.pnas.org/doi/pdf/10.1073/pnas.1918937117


# Environment set up
The following software is required:

*R version 1.9.0 or later with core packages:
  * umap (version 0.2.3.1) 
  * infercnv (version 1.0.2)
  * CellRanger (version 3.0.2)
  * zinbwave (version 1.8.0)
  * SingleR 
  * ImmClassifier
  * SingleCellNet
  * Seurat
  * mclust
  * fastcluster
  * GSVA
  * mgcv
  * lmer4
  * ggplot2
  * data.table
  * tidyr
  * dplyr
  * ggsci
  * vegan
  * boot
  * lmerTest
  * igraph
  * ggraph


# Data availablity
Raw single cell RNA-seq data are available through GEO under accession code GSE211434. Source data are provided with the accompanying manuscript. For each main results figure (Figs 2-6), source data (named: SourceData_Figure#_BriefTitle.csv) is  provided within the source data folder. This input csv file provides the data used in each analysis of the figure. We also provide the output dataset produced by the analysis conducted to generate each subpanel (files in the Outputs subfolder with the suffix “_Output.csv”). Source data has been provided for all main and supplementary figures. Other available data supporting findings are available from the corresponding authors on reasonable request.

These scRNAseq analyses use three main data inputs (detailed below): 
* Tumor microenvironment cell type annotations (composition input: Figure 2/5)
* Cell type gene expression profiles (phenotype input: Figure 4/5)
* Ligand-receptor expression profiles (communication input: Figure 3-5)


## Tumor microenvironment cell type annotation and verification
Stringent quality controls were developed by Griffiths et al. (2021). This was aplied to raw single cell RNA-seq data, to obtain high quality transcriptomic profiles for 424,581 single cells with high-coverage, low mitochondrial content and high confidence of doublet removal. 

Code to conduct these analyses are available here:
https://github.com/U54Bioinformatics/FELINE_project/tree/master/scRNA_cellranger/scripts
https://github.com/U54Bioinformatics/FELINE_project/tree/master/scRNA_doublet

Broad cell types were annotated using singleR, cancer cells were identified by their frequent and pronounced copy numbers amplification using inferCNV. Cell type annotations were verified by cell type specific marker gene expression and UMAP/TSNE analyses. 

Code to conduct these analyses was developed for Grifiths et al. (2021) and are available here:
https://github.com/U54Bioinformatics/FELINE_project/tree/master/scRNA_seurat_celltype_annotation_10x/scripts

Granular immune subtype annotations were obtained using our recently published ImmClassifier machine learning method, which has been validated by flow cytometry comparisons. 

Code to conduct these analyses are available here:
https://github.com/xliu-uth/ImmClassifier


## Cell type gene expression profiling
Cell phenotypes where quantified by their pathway activity in Grifiths et al. (2021), through Gene Set enrichment analysis ssGSEA)
Code to conduct these analyses are available here:
https://github.com/U54Bioinformatics/FELINE_project/tree/master/scRNA_ssGSEA_10x/scripts

Phenotype heterogeneity within cell types was also characterized through UMAP analysis of single cell transcriptional profiles (log(1+CPM)). Myeloid phenotypic heterogeneity was examined in greater detail. Gene expression profiles of myeloid cells provided in the source data for figure 4 are combined with UMAP dimension reduction coordiantes.


## Ligand-receptor (LR) expression and tumor wide communication profiling
Measurements of tumor-wide signaling from diverse non-cancer cell sub-populations and heterogeneous cancer lineages to each receiving cell were obtained using the extended expression-product approach. A worked example of the application of this algorithm to published scRNA data are provided in the folder: "Communication Peripheral blood application Griffiths 2020"

https://github.com/U54Bioinformatics/FELINE_project/tree/master/FELINE_immune_communication/Communication%20Peripheral%20blood%20application%20Griffiths%202020%20

This example uses publically available data from Griffiths et al. (2020) as an input to generate tumor-wide communication scores between cell type populations. Manuscripts source data provides this data in the required input format.A curated ligand-receptor communication database (Ramilowski 2015) defined a set of LR communication pathways (C_jk (x_k,y_jk )) based on known protein-protein interactions. The set of 1444 LR communications measured from FELINE scRNA data are listed in the source data for figure S4.


# Perform analyses within the "Source code" folder
The above listed input data should be accessed via the manuscript source data folder.
Code in the "Source code" folder of this repository performs analyses presented in the manuscript. Code is partitioned into separate scripts to perform analyses relating to each subpanel of the paper's figures. 

Set the working directory, by changing the file path defined at the beginning of each script, to  the location of the source data folder on your machine. Ensure the file path correctly points to the folder containing the  indicated input file (using file.exists(#filepath#)).

Code should not need to be run in the order presented in the figures, as provided input data is sufficient. However, many supplementary figures were generated during the workflow presented in the main figures and are supportive of those findings. Therefore, those analyses remain integrated within the code for the main figures. Code to perform separate supplementary analyses are provided in the "Supplementary Information analyses" folder.





