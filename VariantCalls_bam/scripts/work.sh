
echo "make a copy of bam files of WES data for 24 FELINE patients"
sbatch --array 1-72 cp_files.sh

