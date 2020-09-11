
##################################Install edited version of velocity.R#########################
echo "install edited version of velocyto.R"
# use this version we can export arrow position to text files for Jason
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_velocityR/lib/R/library
/home/jichen/software/BETSY/install/envs/scRNA_velocityR/bin/R
install.packages("callr")
install.packages("devtools")
library(devtools)
install_github("JinfengChen/velocyto.R")

##################################Run commands below to reproduce the results#########################
echo "Input 1: cell ranger results using mRNA reference (with exon-intron structure). For other analysis we used premRNA reference"
#Folder: FEL011046_cell_ranger.mRNA
#the reference version is GRCh38
#betsy_scRNA_10x_QC.sh is the script how I generate these cell ranger folders
mkdir FEL011046_cell_ranger.mRNA
ln -s /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna//scRNA_cellranger/output/FELINE_cellranger_mRNA/COH0*_cell_ranger/cell_ranger/FEL0* ./FEL011046_cell_ranger.mRNA/

echo "Input 2: Jason's UMAP analyses using subset of genes from each Arm"
#Folder: Jason_input/10_umap_visualization/
#2.1, RData Jason save UMAP results: ARM_C_Day180_cancer_cell_UMAP_results.RData
#2.1, gene list Jason used to do UMAP: ARM_A_cancer_cell_UMAP_results.gene.txt
#2.2, UMAP coordinate Jason generated: ARM_A_cancer_cell_UMAP_results.umap.txt
# gene list and UMAP coordinate can be generated from RData using the the following command
bash Jason_input.sh

echo "Step 1: Run velocyto to estimate spliced and unspliced count"
sbatch --array 1-37 scRNA_RNA_velocity_step1.sh

echo "Step 2: merge individual loom file from each patient sample to a single large loom for all FELINE patients"
#patients samples are listed in "samples.list". This can be done for a subset of sample for testing (samples.subset.list). The command below can generate a sample.list file
ls -d FEL011046_cell_ranger.mRNA/*/velocyto/*.loom | awk '{ print "basename "$1" .loom"}' > sample.sh
bash sample.sh > samples.list
#run merging with python "loompy" module
sbatch scRNA_RNA_velocity_step2.sh

echo "Step 3: Run RNA velocity with SeuratWrappers/velocyto.R"
# This will run ArmA, B, C
# Output 1: RNA velocity plot on de novo UMAP coordinate: FELINE_test_ArmA_velocity_denovo_umap.pdf
# Output 2: RNA velocity plot on Jason's UMNAP coordinate: FELINE_test_ArmA_velocity_given_umap.pdf 
# Output 3: cell arrow length and direction predicted by RNA velocity: FELINE_test_ArmA_velocity.arrows.txt
# Output 4: field arrow length and direction predicted by RNA velocity: FELINE_test_ArmA_velocity.garrows.txt
sbatch --array 1-3 scRNA_RNA_velocity_step3.sh


