
#######################commands below are to reproduce the results ###################################
echo "prepare 10x and icell8 ssGSEA score to merge"
ln -s ../scRNA_ssGSEA_10x/output/*.ssGSEA_scores.txt ./
ln -s ../scRNA_ssGSEA_ICELL8/output/*.ssGSEA_scores.txt ./

echo "merge ssGSEA for each cell type: 9 cell type found in 10x data"
sbatch --array 1-9 Normalize_count_by_zinbwave_merge_multiple.sh

echo "convert txt file to RDS file"
sbatch --array 1-9 Matrix2RDS.sh

