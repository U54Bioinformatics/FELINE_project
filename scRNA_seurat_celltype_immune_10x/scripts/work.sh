#######################commands below are to reproduce the results ###################################

echo "prepare scripts and data for refining immunce cells"
ln -s ../scRNA_seurat_integration_10x/output/FEL011046_10x_Seurat_2kgenes_vst_cc.integrated_RPCA.rds FEL011046_10x_Seurat_2kgenes_vst_cc.rds
cp ../scRNA_seurat_integration_10x_umap/input/FEL011046.cell_metadata.txt ./

echo "extract Macrophages, T cells, B cells; did Endothelial cells too"
sbatch --array 1-4 scRNA_Seurat_extract_celltype_obj.sh

echo "recluster and plot markers; did Endothelial cells too"
sbatch --array 1-4 scRNA_Seurat_count2umap.sh

echo "refine immune clusters"
# B cells
## summary cluster annotation; need to manual edit when needed
python scRNA_Seurat_05individual_from_individual_sum_cluster_cell_type.py --meta FEL011046_B_cells_10x.filtered.counts.cell_metadata.UMAPcluster1.txt --cluster_col RNA_snn_res.0.8 --celltype_col Celltype_subtype
## update annotation using cluster annotation above
python scRNA_Seurat_05individual_from_individual_update_cell_type_anno.py --cluster_anno FEL011046_B_cells_10x.filtered.counts.cell_metadata.UMAPcluster1.crosstab_cluster_celltype.cluster_anno.txt --meta_cluster FEL011046_B_cells_10x.filtered.counts.cell_metadata.UMAPcluster1.txt --cluster RNA_snn_res.0.8 --replace 1 --replace_anno Celltype_subtype
## extract short anno for final meta file
python scRNA_Seurat_05individual_from_individual_extract_short_anno.py --input FEL011046_B_cells_10x.filtered.counts.cell_metadata.UMAPcluster1.revised_anno.txt
# T cells
##summary cluster annotation; need to manual edit when needed
python scRNA_Seurat_05individual_from_individual_sum_cluster_cell_type.py --meta FEL011046_T_cells_10x.filtered.counts.cell_metadata.UMAPcluster1.txt  --cluster_col RNA_snn_res.2 --celltype_col Celltype_subtype
## update annotation using cluster annotation above
python scRNA_Seurat_05individual_from_individual_update_cell_type_anno.py --cluster_anno FEL011046_T_cells_10x.filtered.counts.cell_metadata.UMAPcluster1.crosstab_cluster_celltype.cluster_anno.txt --meta_cluster FEL011046_T_cells_10x.filtered.counts.cell_metadata.UMAPcluster1.txt --cluster RNA_snn_res.2 --replace 1 --replace_anno Celltype_subtype
## extract short anno for final meta file
python scRNA_Seurat_05individual_from_individual_extract_short_anno.py --input FEL011046_T_cells_10x.filtered.counts.cell_metadata.UMAPcluster1.revised_anno.txt
# Macrophage
##summary cluster annotation; need to manual edit when needed
python scRNA_Seurat_05individual_from_individual_sum_cluster_cell_type.py --meta FEL011046_Macrophages_10x.filtered.counts.cell_metadata.UMAPcluster1.txt --cluster_col RNA_snn_res.0.8 --celltype_col Celltype_subtype
## update annotation using cluster annotation above
python scRNA_Seurat_05individual_from_individual_update_cell_type_anno.py --cluster_anno FEL011046_Macrophages_10x.filtered.counts.cell_metadata.UMAPcluster1.crosstab_cluster_celltype.cluster_anno.txt --meta_cluster FEL011046_Macrophages_10x.filtered.counts.cell_metadata.UMAPcluster1.txt --cluster RNA_snn_res.0.8 --replace 1 --replace_anno Celltype_subtype
## extract short anno for final meta file
python scRNA_Seurat_05individual_from_individual_extract_short_anno.py --input FEL011046_Macrophages_10x.filtered.counts.cell_metadata.UMAPcluster1.revised_anno.txt
