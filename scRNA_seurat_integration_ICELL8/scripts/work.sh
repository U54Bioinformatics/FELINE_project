
echo "run integration with seurat for 10 icell8 patients: clusters are not good enough to seperate cell types"
# 1. input 1 is FEL0*_10x_Seurat_2kgenes_vst_cc.raw.rds, which include raw cells for each patient
#    input 2 is FEL001010_icell8.cell_type.txt, a raw cell type annotation.
# 2. following the method decribed in https://satijalab.org/seurat/v3.2/integration.html
#    Method used: Reciprocal PCA with few modifications
# 3. output 1: FEL001010_icell8_Seurat_2kgenes_vst_cc.integrated_RPCA.rds.
#    output 2: FEL001010_icell8_Seurat_2kgenes_vst_cc.filtered.rds, merge and filtered cells
sbatch scRNA_Seurat_01cluster_integrate_RPCA_object_icell8.sh

echo "extract counts for seurat to recluster"
# input 1: FEL001010_icell8_Seurat_2kgenes_vst_cc.filtered.rds
# input 2: FEL001010_icell8.cell_metadata.UMAPcluster.txt 
# output 1: FEL001010_icell8_gene_symbols.raw.counts.txt
# output 2: FEL001010_icell8_gene_symbols.CPM.txt
# output 3: FEL001010_icell8.cell_metadata.UMAPcluster.txt 
ln -s FEL001010_icell8_cell_metadata.UMAPcluster.filtered.txt FEL001010_icell8.cell_metadata.UMAPcluster.txt
sbatch scRNA_Seurat_export_count_matrix.sh

echo "do a standard seurat clustering"
# input 1: FEL001010_icell8_gene_symbols.raw.counts.txt
# input 2: FEL001010_icell8.cell_metadata.UMAPcluster.txt
# output : FEL001010_icell8_anno.filtered.counts.seurat_standard.rds. this is the seurat obj used to generate umap
ln -s FEL001010_icell8_gene_symbols.raw.counts.txt FEL001010_icell8_anno.txt
ln -s FEL001010_icell8.cell_metadata.UMAPcluster.txt FEL001010_icell8_anno.metadata.txt
sbatch scRNA_Seurat_count2umap.sh


