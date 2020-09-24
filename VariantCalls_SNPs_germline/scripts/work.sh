
echo "call germline mutations for FEL021 and FEL045, who have higher tumor burden"
# input 1: FELINE5_sample.txt
# input 2: FELINE2_normal_cancer.txt 
# input 3: bam files
mkdir FELINE2
ln -s /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_bam/output/FELINE_bam/FEL021_* FELINE2/
ln -s /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_bam/output/FELINE_bam/FEL045_* FELINE2/
# make sure there is one line with "FELINE2" in file "FELINE_patients.4.list"
# in this case "FELINE2" is the second line, so we use SLURM with "--array 2".
sbatch --array 2 betsy_WGS_germline_variants_FELINE.sh

echo "cal germline mutations for all 24 patients"
ln -s /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_bam/output/FELINE_bam ./
sbatch --array 1 betsy_WGS_germline_variants_FELINE.sh

