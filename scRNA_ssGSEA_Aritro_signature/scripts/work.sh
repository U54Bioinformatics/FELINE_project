
echo "convert signatures in RData to gmt format" 
cat AN_GeneSets.R | R --slave

echo "run ssGSEA for 10x zinbwave normalized count" 
ln -s ../scRNA_ssGSEA_10x/output/FEL011046_Cancer_cells_10x.zinbwave.normalized.txt ./
ln -s ../scRNA_ssGSEA_10x/output/FEL011046_Normal_epithelial_cells_10x.zinbwave.normalized.txt ./
sbatch --array 1-2 betsy_scRNA_ssGSEA.sh

echo "run ssGSEA for 10x zinbwave normalized count"
ln -s ../scRNA_ssGSEA_ICELL8/output/FEL001010_Cancer_cells_icell8.zinbwave.normalized.txt ./
ln -s ../scRNA_ssGSEA_ICELL8/output/FEL001010_Normal_epithelial_cells_icell8.zinbwave.normalized.txt ./
sbatch --array 3-4 betsy_scRNA_ssGSEA.sh

echo "merge 10x and icell8"
sbatch --array 1-2 Normalize_count_by_zinbwave_merge_multiple.sh


