##################Run command below to reproduce results##################

##################Process 10x data########################################
echo "prepare seurat obj"
ln -s ../scRNA_seurat_obj/FEL011046_10x_Seurat_2kgenes_vst_cc.rds FEL011046_10x.rds

echo "meta file: The meta file is the original one I used when I do zinbwave normalization. This file has 174066 cells, while the final cell meta has 173160 cells (../scRNA_meta/output/FEL001046_scRNA.metadata.clinical.txt). Some cells are removed later during the analysis due to quality control but no cell was added after that"
#FEL011046_10x.metadata.txt

echo "Prepare count and meta for each cell type. Jason suggests to do zinbwave normalization on each cell type"
# extract Cell.ID for each cell type
# output is FEL011046_*_10x.Cell_ID.txt for each cell type
python Prepare_celltype_cell_list.py --meta FEL011046_10x.metadata.txt --output FEL011046 --platform 10x
# extract meta file from master meta file: FEL011046_10x.metadata.txt
# output is FEL011046_*_10x.metadata.txt for each cell type
sbatch --array 1-9 Extract_cell_data.sh
# extract count file (for cell type except cancer cells) from master rds object: FEL011046_10x.rds
# output is 
sbatch --array 1-8 scRNA_Seurat_export_count_matrix.celltype.sh
# extract count for cancer cells
sbatch --array 1-35 scRNA_Seurat_export_count_matrix.split.sh
#uncomment the fellowing command in "scRNA_Seurat_export_count_matrix.split.sh"
#python scRNA_Seurat_export_count_matrix.split.merge.py --list $prefix\.patients.list --prefix $prefix
sbatch scRNA_Seurat_export_count_matrix.split.sh
mv FEL011046_10x_gene_symbols.raw.counts.txt FEL011046_Cancer_cells_10x.count.txt
mv FEL011046_10x_gene_symbols.CPM.txt FEL011046_Cancer_cells_10x.CPM.txt

echo "run zinbwave normalization"
# output are "FEL011046_*_10x.zinbwave.normalized.txt", "FEL011046_*_10x.zinbwave.raw.txt" and "FEL011046_*_10x.zinbwave.residuals.txt"
# T cells and B cells, macrophage
sbatch --array 1-3 Normalize_count_by_zinbwave.01_prepare.sh
#Bcell
sbatch --array 2 Normalize_count_by_zinbwave.02_normalization_single_file.sh
#Tcell
sbatch --array 1-2 Normalize_count_by_zinbwave.02_normalization.sh
#macrophage
sbatch --array 1-6 Normalize_count_by_zinbwave.02_normalization.sh
sbatch --array 1-3 Normalize_count_by_zinbwave.03_merge_output.sh

# Adipocytes and Pericytes
sbatch --array 4-5 Normalize_count_by_zinbwave.01_prepare.sh
sbatch --array 4-5 Normalize_count_by_zinbwave.02_normalization_single_file.sh
sbatch --array 4-5 Normalize_count_by_zinbwave.03_merge_output.sh

# Endothelial cell
sbatch --array 6 Normalize_count_by_zinbwave.01_prepare.sh
sbatch --array 1-4 Normalize_count_by_zinbwave.02_normalization.sh
sbatch --array 6 Normalize_count_by_zinbwave.03_merge_output.sh

# fibroblast
sbatch --array 7 Normalize_count_by_zinbwave.01_prepare.sh
sbatch --array 1-12 Normalize_count_by_zinbwave.02_normalization.sh
sbatch --array 7 Normalize_count_by_zinbwave.03_merge_output.sh

# normal epithelial cells
sbatch --array 8 Normalize_count_by_zinbwave.01_prepare.sh
sbatch --array 1-8 Normalize_count_by_zinbwave.02_normalization.sh
sbatch --array 8 Normalize_count_by_zinbwave.03_merge_output.sh

# cancer cells
sbatch --array 9 Normalize_count_by_zinbwave.01_prepare.sh
sbatch --array 1-56 Normalize_count_by_zinbwave.02_normalization.sh
sbatch --array 9 Normalize_count_by_zinbwave.03_merge_output.sh



