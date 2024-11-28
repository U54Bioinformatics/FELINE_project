# Blocking cancer-fibroblast mutualism inhibits proliferation of endocrine therapy resistant breast cancer

This repository contains code accompanying the manuscript: "Blocking cancer-fibroblast mutualism inhibits proliferation of endocrine therapy resistant breast cancer". It provides code used to analyze the breast tumor cell heterogeneity, identify compensatory growth factor signaling in cancer cells resistant to therapy and decipher tumor-wide communications from  diverse populations of cancer and non-cancer cells within the tumor microenvironment (TME). Serial ER+ breast tumor biopsies from of post-menopausal women on the FELINE clinical trial (clinicaltrials.gov # NCT02712723) were profiled using single cell RNA-seq (scRNAseq). Tumor samples from two cohorts of patients were independently examined, with the first discovery cohort used to identify tumor features promoting treatment resistance. The second cohort was used to validate discoveries and critical predictions. Patient derived predictions about the nature and source of compensatory growth factor signaling were then tested across multiple in vitro cell line model systems. Coculutre experiments were used to verify the role of non-cancer cell types (fibroblasts) in supplying growth factor when coculutred with cancer cells or when under drug pressure. Targeted inhibitors of fibroblast-cancer communications were shown to block fibroblasts from facilitating cancer growth and to effectively control the growth of endocrine +/- CDKi sensitive and resistant cancer cell lines.

Related publications (Please consult and cite accordingly when using this repository):

Griffiths, J. I., Chen, J., Cosgrove, P. A., O’Dea, A., Sharma, P., Ma, C., ... & Bild, A. H. (2021). Serial single-cell genomics reveals convergent subclonal evolution of resistance as patients with early-stage breast cancer progress on endocrine plus CDK4/6 therapy. Nature cancer, 2(6), 658-671.https://www.nature.com/articles/s43018-021-00215-7

Griffiths, J. I., Wallet, P., Pflieger, L. T., Stenehjem, D., Liu, X., Cosgrove, P. A., ... & Bild, A. H. (2020). Circulating immune cell phenotype dynamics reflect the strength of tumor–immune cell interactions in patients during immunotherapy. Proceedings of the National Academy of Sciences, 117(27), 16072-16082. https://www.pnas.org/doi/pdf/10.1073/pnas.1918937117

Griffiths, J. I., Cosgrove, P. A., Castaneda, E. M., Nath, A., Chen, J., Adler, F. R., ... & Bild, A. H. (2022). Cancer cells communicate with macrophages to prevent T cell activation during development of cell cycle therapy resistance. bioRxiv, 2022-09.


# Environment set up
The following software is required:

*R version 1.9.0 or later with core packages:
  * zinbwave (version 1.8.0)
  * infercnv (version 1.0.2)
  * CellRanger (version 3.0.2)
  * SingleR 
  * ImmClassifier
  * SingleCellNet
  * umap (version 0.2.3.1) 
  * Seurat
  * mclust
  * fastcluster
  * Rdimtools
  * Rfast
  * GSVA
  * mgcv
  * lmer4
  * lmerTest
  * merTools
  * effects
  * ggplot2
  * data.table
  * tidyr
  * dplyr
  * ggsci
  * viridis
  * vegan
  * boot
  * igraph
  * ggraph
  * ggdendro
  * dendextend
  * ider


# Data availablity
Raw single cell RNA-seq data are available through GEO under accession code GSE211434. Source data are provided with the accompanying manuscript. For each main results figure (Figs 2-6), source data (named: SourceData_Figure#_BriefTitle.csv) is provided within the source data folder. This input csv file provides the data used in each analysis of the figure. We also provide image data for representative images presented in the main text and the experimental replicate images supporting or findings.Other available data supporting findings are available from the corresponding authors on reasonable request.

These scRNAseq analyses use three main data inputs (detailed below): 
* Tumor microenvironment cell type annotations to distinguish cancer cells (Figure 2) and their communication with other cancer and non-cancer cell types (Figure 3A/B and 4A/B)
* Cell type gene expression profiles and gene set enrichment (ssGSEA) scores (Figure 2A and 3C)
* Ligand-receptor expression profiles (cell type communication: Figure 2-4)

In vitro validation experiments provided three main data type:
* mRNA expression (Figure 2A,3E and 4C/D)
* western blot protein quantification (Figure 2C)
* co-cultured cancer and non-cancer population growth under treatment (Figure 5-6)

## Tumor microenvironment cell type annotation and verification
Rigorous quality control (QC) measures, implemented by Griffiths et al. (2021), were applied to raw scRNAseq data, to obtain high quality transcriptomic profiles for 424,581 single cells. These scRNAseq profiles have high-coverage, low mitochondrial content and high confidence of doublet removal. 

Code to conduct these QC analyses are provided here:
https://github.com/U54Bioinformatics/FELINE_project/tree/master/scRNA_cellranger/scripts
https://github.com/U54Bioinformatics/FELINE_project/tree/master/scRNA_doublet

Broad cell types were annotated using singleR. Cancer cells were identified by their frequent and pronounced copy numbers amplification using inferCNV. Cell type annotations were verified by cell type specific marker gene expression and UMAP/TSNE analyses. 

Code to conduct these cell type analyses was developed for Griffiths et al. (2021) and are available here:
https://github.com/U54Bioinformatics/FELINE_project/tree/master/scRNA_seurat_celltype_annotation_10x/scripts

Granular immune subtype annotations were obtained using the ImmClassifier machine learning method.

Code to conduct these analyses are available here:
https://github.com/xliu-uth/ImmClassifier

Cell type and subtype annotations for each single cell across the temporal patient tumor samples from the discovery and validation data sets are provided in source data of Griffiths et al. (2022). Cells within the dataset annotated as fibroblasts or cancer cells are listed in source data for figure 2A and 3C.


## Cell type gene expression profiling
Cell phenotypes where quantified by their pathway activity in Griffiths et al. (2021), through Gene Set enrichment analysis (ssGSEA).
Code to conduct these analyses are available here:
https://github.com/U54Bioinformatics/FELINE_project/tree/master/scRNA_ssGSEA_10x/scripts

Cancer ERBB growth factor pathway activation during treatment was a key focused of this study (See source data for Figure 2A). Phenotype heterogeneity within cell types was also characterized through UMAP analysis of scRNAseq transcriptional profiles (log(1+CPM)). ERBB ligand gene expression profiles of fibroblast cells are combined with UMAP dimension reduction coordinates (provided source data for figure 3C).


## Ligand-receptor (LR) expression and tumor wide communication profiling
Measurements of tumor-wide signaling from diverse non-cancer cell sub-populations and heterogeneous cancer lineages to each receiving cell were obtained using the extended expression-product approach (see Armingol et al. 2021) . A worked example of the application of this algorithm to published scRNA data are provided here: 
https://github.com/U54Bioinformatics/FELINE_project/tree/master/FELINE_immune_communication/Communication%20Peripheral%20blood%20application%20Griffiths%202020%20

This example uses publicly available data from Griffiths et al. (2020) as an input to generate tumor-wide communication scores between cell type populations. Manuscripts source data provides this data in the required input format.A curated Ligand-receptor communication database (Ramilowski 2015) defined a set of LR communication pathways based on known protein-protein interactions. The set of LR communications measured from FELINE scRNA data are listed in the source data for figures 3A and 4A.


## mRNA expression of cancer and fibrobalst cells
Quantitative real time-PCR analysis measured treatment induced gene expression changes under treatment controlled experimental conditions (Source data for Figure 2A,3E and 4C/D). Replicate mRNA measurements are normalized relative to treatment controls.


## Western blot protein quantification 
Immunoblotting was applied to assess  endocrine (fulvestrant) induced ERBB protein phosphorylation. Replicate experiments were performed and western blot bands were quantified using imageJ (see Source data for Figure 2C).


## Co-cultured cancer and non-cancer population growth rate 
To assess the ability of fibroblasts to facilitate cancer growth, we co-cultured these cell types (fluorescently labelled) in 3D spheroids and measured log linear population growth over time (relative growth rate:RGR) using repeated fluorescence imaging (Figure 5-6)


# Perform analyses within the "Source code" folder
The above mentioned input data should be accessed via the manuscript source data folder (with Raw data at GEO GSE211434). 
Code in the "Source code" folder of this repository performs analyses presented in the manuscript. Code is partitioned into separate scripts which each perform an analysis relating to a specific subpanel within a manuscript figure. Code names begin with a string that indicates the figure panel.

Edit the file path defined at the beginning of each script ("SourceDataLoc") to match the location of the source data folder on your machine. Ensure the file path correctly points to the folder containing the indicated input file (using file.exists(#filepath#)). Packages required are listed at the top of each script and need installation prior to running. Code to save output is currently commented out (users should indicate appropriate output file locations).

Scripts should not need to be run in the order presented in the figures, as the provided input data is sufficient. However, supplementary figures were generated during the workflow presented in the manuscript and are supportive of those findings. Therefore, those analyses remain integrated within the code for the main figures. 





