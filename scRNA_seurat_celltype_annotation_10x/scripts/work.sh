

echo "add infercnv and scrublet (doublet) to meta"
# input 1: FEL011046_10x_cell_metadata.UMAPcluster.integrated_RPCA.txt, cell meta from seurat intergration
# input 2: FEL011046_cell_metadata.infercnv_CNA.txt, infercnv call CNA or noCNA
# input 3: FEL011046.scrublet_out.csv, doublet prediction
# output: FEL011046_10x_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.txt
cp ../scRNA_seurat_integration_10x/output/FEL011046_10x_cell_metadata.UMAPcluster.integrated_RPCA.txt .
cp ../scRNA_doublet/output/FEL011046.scrublet_out.csv ./
cp ../scRNA_inferCNV_heatmap/output/FEL011046_cell_metadata.infercnv_CNA.txt ./
python scRNA_Seurat_05individual_from_individual_add_infercnv_CNA_anno.py --meta_cluster FEL011046_10x_cell_metadata.UMAPcluster.integrated_RPCA.txt --infercnv FEL011046_cell_metadata.infercnv_CNA.txt
python scRNA_Seurat_05individual_from_individual_add_scrublet_anno.py --meta_cluster FEL011046_10x_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.txt --scrublet FEL011046.scrublet_out.csv

echo "Annotate cell type for each seurat cluster based on singleR and marker gene expression"
# config running parameters
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
sample=FEL011046
platform=10x
prefix=$sample\_$platform
# generate a seurat cluster vs. sample matrix for QC
# output file is FEL011046_10x_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.sample_matrix.txt 
/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_sample_cluster_matrix.R $prefix
# generate a seurat cluster vs. singleR annotation matrix for QC
# output file is FEL011046_10x_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.SingleR_cluster_table_Blueprint.txt 
/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_05individual_from_individual_singleR_cluster_matrix.R $prefix
# summarize singleR annotation for each cluster and select the annotation that has maximum cells in each cluster
# output file is FEL011046_10x_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.SingleR_cluster_table_Blueprint.cluster_anno.txt 
# This output file is where we can use marker gene expression to validate annotation of each cluster and make changes to this file if the annotation is incorrect.
# For 10x data, I annotated patients several times. So I used previous annotation that based on singleR and marker gene expression to summerize cluster annotation.
python scRNA_Seurat_05individual_from_individual_sum_cell_type.py --cluster_anno $prefix\_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.SingleR_cluster_table_Blueprint.txt --sample_cluster $prefix\_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.sample_matrix.txt
# use the finalized *.cluster_anno.txt file above to update the cell meta file: FEL011046_10x_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.txt
# output file is FEL011046_10x_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.revised_anno.txt
# I already annotated P11-43 in previous run, so for final results I only updated P44-46 when the data comes. 
python scRNA_Seurat_05individual_from_individual_add_cell_type_anno_update044046.py --cluster_anno $prefix\_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.SingleR_cluster_table_Blueprint.cluster_anno.txt --meta_cluster $prefix\_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.txt

# the annotation from this analysis is good enough for broad range cell type annotation. To further subclassify B cells, T cells, Macrophage, and Endothelial cells, I recluster these cells in folder "scRNA_seurat_celltype_immune_10x" and these annotation were added to "scRNA_meta/output/FEL001046_scRNA.metadata.clinical.txt"

