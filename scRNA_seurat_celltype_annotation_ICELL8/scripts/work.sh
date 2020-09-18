#######################commands below are to reproduce the results ###################################


echo "Annotate cell type for each seurat cluster based on singleR and marker gene expression"
# output is FEL001010_icell8_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.revised_anno.txt
cp ../scRNA_seurat_integration_ICELL8/output/FEL001010_icell8_anno.filtered.counts.cell_metadata.UMAPcluster1.txt FEL001010_icell8_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.txt
# config running parameters
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
sample=FEL001010
platform=icell8
prefix=$sample\_$platform
# generate a seurat cluster vs. sample matrix for QC
# output file is FEL001010_icell8_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.sample_matrix.txt
/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_sample_cluster_matrix.R $prefix
# generate a seurat cluster vs. singleR blueprint cell annotation (column Encode_main_type in input file) matrix for annotation
# output file is FEL001010_icell8_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.SingleR_cluster_table_Blueprint.txt 
/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_singleR_cluster_matrix.R $prefix
# summarize singleR annotation for each cluster and select the annotation that has maximum cells in each cluster
# output file is FEL001010_icell8_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.SingleR_cluster_table_Blueprint.cluster_anno.txt
# Cluster	Max	CellType
# 0	136	CD8+ T cells
# 1	118	Macrophages
# 2	247	Epithelial cells
# This output file is where we can use marker gene expression to validate annotation of each cluster and make changes to this file if the annotation is incorrect. 
python scRNA_Seurat_05individual_from_individual_sum_cell_type.py --cluster_anno $prefix\_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.SingleR_cluster_table_Blueprint.txt --sample_cluster $prefix\_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.sample_matrix.txt
# use the finalized *.cluster_anno.txt file above to update the cell meta file: FEL001010_icell8_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.txt
# output file is FEL001010_icell8_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.revised_anno.txt 
python scRNA_Seurat_05individual_from_individual_add_cell_type_anno.py --cluster_anno $prefix\_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.SingleR_cluster_table_Blueprint.cluster_anno.txt --meta_cluster $prefix\_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.txt
# when I rerun this analysis I notice there is slight difference in cluster numbers between this run my original run. So the reproduced annotation might be slightly difference from the original (final) one. I put both original run and reproduce run in folders: FEL001010_original_cluster_anno and FEL001010_reproduce_cluster_anno. The results "FEL001010_original_cluster_anno" was used in following analysis.


echo "Update cell type annotation with seurat cluster based annotation: epithelial cells, fibroblast cells and endothelial cells are identifed; Immune cells are not subclassified yet;"
# input 1 is FEL001010_icell8.cell_metadata.UMAPcluster.txt
# input 2 is FEL001010_icell8_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.revised_anno.txt
# output is FEL001010_icell8.cell_metadata.UMAPcluster.revised_anno.txt
# in the output file, column "Encode_main_cluster" has updated cell type annotation
python scRNA_Seurat_reformat_icell8_annotation.py --meta FEL001010_icell8.cell_metadata.UMAPcluster.txt --cell_anno FEL001010_icell8_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.revised_anno.txt

echo "extract immune cells and fine annotation immune cells"
# this is done in a different folder: scRNA_seurat_celltype_immune_ICELL8
cp ../scRNA_seurat_celltype_immune_ICELL8/output/FEL001010_icell8.immune.metadata.revised_anno.txt ./

echo "Finalize cell type annotaiton with immune cell annotation and infercnv calls"
# input 1: updated cell type annotation from above step
# input 2: FEL001010_cell_metadata.infercnv_tumor_normal.txt, infercnv calls (CNA and noCNA indicating CNA status)
# input 3: FEL001010_icell8.immune.metadata.revised_anno.txt, immune cell annotation with B cells, T cells and Macrophage annotated
# output file is FEL001010_icell8.cell_metadata.UMAPcluster.revised_anno.finalized_short.txt, which is the final annotation for ICELL8 data: 3736 cells in this file = 3484 final cells + 252 low-quality cells
# infercnv calls is used to 1, distinguish cancer cells and normal epithelial cells; 2, to assign immune cells aor stromal cells to low-quality if they have CNA.
python scRNA_Seurat_finalize_icell8_annotation.py --meta FEL001010_icell8.cell_metadata.UMAPcluster.revised_anno.txt --infercnv_anno FEL001010_cell_metadata.infercnv_tumor_normal.txt --immune_anno FEL001010_icell8.immune.metadata.revised_anno.txt


