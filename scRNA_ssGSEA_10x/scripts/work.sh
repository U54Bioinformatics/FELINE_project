
echo "prepare zinbwave normalized counts"
ln -s ../scRNA_count_zinbwave_10x/output/FEL011046_*_10x.zinbwave.normalized.txt ./

echo "run ssGSEA with BETSY"
sbatch --array 1-9 betsy_scRNA_ssGSEA.sh

