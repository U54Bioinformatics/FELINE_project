
#######################commands below are to reproduce the results ###################################
echo "prepare 10x and icell8 zinbwave count to merge"
ln -s ../scRNA_count_zinbwave_10x/output/FEL011046_*.zinbwave.normalized.txt ./
ln -s ../scRNA_count_zinbwave_ICELL8/output/FEL001010_*.zinbwave.normalized.txt ./

echo "merge zinbwave count for each cell type: 9 cell type found in 10x data"
sbatch --array 1-9 Normalize_count_by_zinbwave_merge_multiple.sh

echo "convert txt file to RDS file"
sbatch --array 1-9 Matrix2RDS.sh

