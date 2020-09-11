
echo "prepare zinbwave normalized counts"
ln -s ../scRNA_count_zinbwave_ICELL8/output/FEL001010_*_icell8.zinbwave.normalized.txt ./

echo "run ssGSEA with BETSY"
sbatch --array 1-10 betsy_scRNA_ssGSEA.sh

