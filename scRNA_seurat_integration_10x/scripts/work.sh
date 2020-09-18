

echo "run integration with seurat for 35 10x patients"
# 1. input 1 is FEL0*_10x_Seurat_2kgenes_vst_cc.raw.rds, which include raw cells for each patient
#    input 2 is FEL011046_10x.cell_type.txt, a raw cell type generated based singleR and marker gene expression. As#    we processed patients from multiple runs, the results are combined from previous runs. This result is not nece#    ssary accurate. We refined cell type annotation after clustering. 
# 2. following the method decribed in https://satijalab.org/seurat/v3.2/integration.html
#    Method used: Reciprocal PCA with few modifications
# 3. output is FEL011046_10x_Seurat_2kgenes_vst_cc.integrated_RPCA.rds. This file is used for final analyses and fi     gures. Only need to update cell meta file when necessary.  
sbatch scRNA_Seurat_01cluster_integrate_RPCA_object.sh

echo "Final seurat obj with 35 patients with 10x scRNA data"
ln -s FEL011046_10x_Seurat_2kgenes_vst_cc.integrated_RPCA.rds FEL011046_10x_Seurat_2kgenes_vst_cc.rds
