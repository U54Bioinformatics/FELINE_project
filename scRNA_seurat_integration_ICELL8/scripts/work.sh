#######################commands below are to reproduce the results ###################################

# For icell8, the Seurat intergration pipeline does not work as expected. A standard merge and clustering was used in the final analysis.
echo "merge and filter 10 icell8 patients"
# input 1 is FEL0*_10x_Seurat_2kgenes_vst_cc.raw.rds, which include raw cells for each patient
# input 2 is FEL001010_icell8.cell_type.txt, a raw cell type annotation collected from singleR
# output 1 is FEL001010_icell8_cell_metadata.UMAPcluster.filtered.txt, cell meta
# output 2 is FEL001010_icell8_Seurat_2kgenes_vst_cc.filtered.rds, seurat obj 
sbatch scRNA_Seurat_01cluster_integrate_RPCA_object_icell8.sh

echo "extract counts for seurat to recluster"
# input 1: FEL001010_icell8_Seurat_2kgenes_vst_cc.filtered.rds
# output 1: FEL001010_icell8_gene_symbols.counts.txt
# output 2: FEL001010_icell8_gene_symbols.CPM.txt
# output 3: FEL001010_icell8.cell_metadata.UMAPcluster.txt 
ln -s FEL001010_icell8_Seurat_2kgenes_vst_cc.filtered.rds FEL001010_icell8_Seurat_2kgenes_vst_cc.rds
sbatch scRNA_Seurat_export_count_matrix.sh

echo "do a standard seurat clustering"
# input 1: FEL001010_icell8_gene_symbols.counts.txt
# input 2: FEL001010_icell8.cell_metadata.UMAPcluster.txt
# output : FEL001010_icell8_anno.filtered.counts.seurat_standard.rds. this is the seurat obj used to generate umap
ln -s FEL001010_icell8_gene_symbols.counts.txt FEL001010_icell8_anno.txt
ln -s FEL001010_icell8.cell_metadata.UMAPcluster.txt FEL001010_icell8_anno.metadata.txt
sbatch scRNA_Seurat_count2umap.sh


