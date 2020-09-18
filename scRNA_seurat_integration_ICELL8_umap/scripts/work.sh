
#######################commands below are to reproduce the results ###################################
echo "prepare input"
# seurat obj
ln -s ../scRNA_seurat_integration_ICELL8/output/FEL001010_icell8_anno.filtered.counts.seurat_standard.rds ./
ln -s FEL001010_icell8_anno.filtered.counts.seurat_standard.rds FEL001010_integrated_RPCA_icell8_Seurat_2kgenes_vst_cc.rds
# cell type annotaiton
cp ../scRNA_meta/output/FEL001046_scRNA.metadata.clinical.txt FEL001046_scRNA.metadata.clinical.txt
python FEL001010.cell_metadata.py --input FEL001046_scRNA.metadata.clinical.txt
ln -s FEL001010.cell_metadata.txt FEL001010_integrated_RPCA_icell8.broadclass.metadata.txt

echo "plot umap cell type annotaiton and marker gene expression"
# Supplementary Figure 1a: plot umap for all cells from 10 patients
# output is FEL001010_integrated_RPCA_icell8_Seurat_UMAP_Celltype/FEL001010_integrated_RPCA_icell8_Celltype.ggplot.broadclass.png
sbatch scRNA_Seurat_05individual_from_individual_featureplots_pdf.sh
# Supplementary Figure 1c: plot panmarker gene expression
# output is FEL001010_integrated_RPCA_icell8_Seurat_marker_genes/FEL001010_integrated_RPCA_icell8_Cell_Type_panmarker_feature.tsne.png
sbatch scRNA_Seurat_05individual_from_individual_plot_marker_genes.sh
# Supplementary Figure 1d: plot barplot for cell type in each patient
# output is FEL001010_scRNA.celltype.barplot.pdf 
bash FEL001010_scRNA.celltype_barplot.sh

