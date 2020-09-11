
#######################commands below are to reproduce the results ###################################
echo "prepare input"
# seurat obj
ln -s ../scRNA_seurat_integration_10x/output/FEL011046_10x_Seurat_2kgenes_vst_cc.integrated_RPCA.rds ./
ln -s FEL011046_10x_Seurat_2kgenes_vst_cc.integrated_RPCA.rds FEL011046_integrated_RPCA_10x_Seurat_2kgenes_vst_cc.rds
# cell type annotaiton
cp ../scRNA_meta/output/FEL001046_scRNA.metadata.clinical.txt FEL001046_scRNA.metadata.clinical.txt
python FEL011046.cell_metadata.py --input FEL001046_scRNA.metadata.clinical.txt
ln -s FEL011046.cell_metadata.txt FEL011046_integrated_RPCA_10x.broadclass.metadata.txt

echo "plot umap cell type annotaiton and marker gene expression"
# Figure 1b: plot umap for all cells from 35 patients
# output is FEL011046_integrated_RPCA_10x_Seurat_UMAP_Celltype/FEL011046_integrated_RPCA_10x_Celltype.ggplot.broadclass.png
sbatch scRNA_Seurat_05individual_from_individual_featureplots_pdf.sh
# Figure 1c: plot panmarker gene expression
# output is FEL011046_integrated_RPCA_10x_Seurat_marker_genes/FEL011046_integrated_RPCA_10x_Cell_Type_panmarker_feature.tsne.png 
sbatch scRNA_Seurat_05individual_from_individual_plot_marker_genes.sh
# Figure 1d: plot barplot for cell type in each patient
# output is FEL001046_scRNA.celltype.barplot.pdf 
bash FEL001046_scRNA.celltype_barplot.sh
