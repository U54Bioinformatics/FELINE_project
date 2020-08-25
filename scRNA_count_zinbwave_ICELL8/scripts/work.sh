##################Run command below to reproduce results##################

##################Process ICELL8 data########################################
echo "master meta and count files"
#FEL001010_icell8.metadata.txt
#FEL001010_icell8.count.txt

echo "Prepare count and meta for each cell type. Jason suggests to do zinbwave normalization on each cell type"
# extract Cell.ID for each cell type
python Prepare_celltype_cell_list.py --meta FEL001010_icell8.metadata.txt --output FEL001010 --platform icell8
# extract meta and count files from master files: FEL001010_icell8.metadata.txt and FEL001010_icell8.count.txt
sbatch --array 1-10 Extract_cell_data.sh


echo "run normalization"
# Run for each cell type by array job using FEL001010_icell8.celltype.txt
# For icell8, all cell types have fewer cells, we can run all by array job
# output are "FEL001010_*_ICELL8.zinbwave.normalized.txt", "FEL001010_*_ICELL8.zinbwave.raw.txt" and "FEL001010_*_ICELL8.zinbwave.residuals.txt"
sbatch --array 1-10 Normalize_count_by_zinbwave.01_prepare.sh
sbatch --array 1-10 Normalize_count_by_zinbwave.02_normalization_single_file.sh
sbatch --array 1-10 Normalize_count_by_zinbwave.03_merge_output.sh




