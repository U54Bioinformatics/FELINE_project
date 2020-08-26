##################Run command below to reproduce results##################

##################Process 10x data########################################
echo "prepare seurat obj"
ln -s ../scRNA_seurat_obj/FEL011046_10x_Seurat_2kgenes_vst_cc.rds ./

echo "extract raw count and CPM for each patient"
sbatch --array 1-2 scRNA_Seurat_export_count_matrix.split.sh

echo "merge data from all patients into a single file"
#this is optional. only useful when a merged raw count or CPM is needed
#uncomment the fellowing command in "scRNA_Seurat_export_count_matrix.split.sh"
#python scRNA_Seurat_export_count_matrix.split.merge.py --list $prefix\.patients.list --prefix $prefix
sbatch scRNA_Seurat_export_count_matrix.split.sh

##################Process ICELL8 data########################################

