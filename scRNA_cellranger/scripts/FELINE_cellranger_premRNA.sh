#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=FELINE_cellranger_premRNA.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./


start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "CPU: $CPU"
echo "N: $N"

echo "COH039: FEL011"
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH039/FEL011.premRNA/FEL011_count03/ ./FELINE_cellranger_premRNA/COH039_FEL011_cell_ranger/cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH039/FEL011.premRNA/FEL011_meta/ ./FELINE_cellranger_premRNA/COH039_FEL011_cell_ranger/qc
#cp /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH039/FEL011.premRNA/FEL011_count01.noaggr.txt ./FELINE_cellranger_premRNA/COH039_FEL011_cell_ranger/read_counts.txt

echo "COH044: FEL012,13,14,15,16"
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL012_count03/ ./FELINE_cellranger_premRNA/COH044_FEL012_cell_ranger/cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL012_meta/ ./FELINE_cellranger_premRNA/COH044_FEL012_cell_ranger/qc
#cp /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL012_count01.noaggr.txt ./FELINE_cellranger_premRNA/COH044_FEL012_cell_ranger/read_counts.txt

#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL013_count03/ ./FELINE_cellranger_premRNA/COH044_FEL013_cell_ranger/cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL013_meta/ ./FELINE_cellranger_premRNA/COH044_FEL013_cell_ranger/qc
#cp /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL013_count01.noaggr.txt ./FELINE_cellranger_premRNA/COH044_FEL013_cell_ranger/read_counts.txt

#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL014_count03/ ./FELINE_cellranger_premRNA/COH044_FEL014_cell_ranger/cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL014_meta/ ./FELINE_cellranger_premRNA/COH044_FEL014_cell_ranger/qc
#cp /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL014_count01.noaggr.txt ./FELINE_cellranger_premRNA/COH044_FEL014_cell_ranger/read_counts.txt

#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL015_count03/ ./FELINE_cellranger_premRNA/COH044_FEL015_cell_ranger/cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL015_meta/ ./FELINE_cellranger_premRNA/COH044_FEL015_cell_ranger/qc
#cp /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL015_count01.noaggr.txt ./FELINE_cellranger_premRNA/COH044_FEL015_cell_ranger/read_counts.txt

#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL016_count03/ ./FELINE_cellranger_premRNA/COH044_FEL016_cell_ranger/cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL016_meta/ ./FELINE_cellranger_premRNA/COH044_FEL016_cell_ranger/qc
#cp /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.premRNA/FEL016_count01.noaggr.txt ./FELINE_cellranger_premRNA/COH044_FEL016_cell_ranger/read_counts.txt

echo "COH049: FEL017,19,20,21,22"
#cp -R /home/jichen/abild/U54_data_folder/COH049/feline_p1722_premrna.noaggr.counts.txt.gz ./FELINE_cellranger_premRNA/COH049_cell_ranger/read_counts.txt.gz
#cp -R /home/jichen/abild/U54_data_folder/COH049/feline_p1722_premrna.metadata ./FELINE_cellranger_premRNA/COH049_cell_ranger/qc
#cp -R /home/jichen/abild/U54_data_folder/COH049/feline_p1722_premrna.10x.counts ./FELINE_cellranger_premRNA/COH049_cell_ranger/cell_ranger

echo "COH047: FEL023,34"
#cp -R /home/jichen/abild/U54_data_folder/COH047/feline_p2324_premrna.noaggr.counts.txt ./FELINE_cellranger_premRNA/COH047_cell_ranger/read_counts.txt
#cp -R /home/jichen/abild/U54_data_folder/COH047/feline_p2324_premrna.metadata/ ./FELINE_cellranger_premRNA/COH047_cell_ranger/qc
#cp -R /home/jichen/abild/U54_data_folder/COH047/feline_p2324_premrna.10x.counts/ ./FELINE_cellranger_premRNA/COH047_cell_ranger/cell_ranger

echo "COH055: FEL025,26,27,28"
cp -R /home/jichen/abild/U54_data_folder/COH055/preprocess_premRNA/ ./FELINE_cellranger_premRNA/COH055_cell_ranger

echo "COH057: FEL029,30,31,32,33,34,35,36,37,38,39,40"
#cp -R /home/jichen/abild/U54_data_folder/COH057/feline_p29_40_premRNA/ ./FELINE_cellranger_premRNA/COH057_cell_ranger

echo "COH064: FEL041,42,43"
#mv /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_10x_scRNA_GATK/FEL011046_cell_ranger_premRNA/COH064_cell_ranger/ ./FELINE_cellranger_premRNA/cell_ranger
#cp -R /home/jichen/abild/U54_data_folder/COH064/Preprocess\ FELINE\ 41-43\ premRNA\ -\ Results/qc/ FELINE_cellranger_premRNA/COH064_cell_ranger/
#cp -R /home/jichen/abild/U54_data_folder/COH064/Preprocess\ FELINE\ 41-43\ premRNA\ -\ Results/read_counts.txt.gz FELINE_cellranger_premRNA/COH064_cell_ranger/

echo "COH064: FEL044,45,46"
#mv /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_10x_scRNA_GATK/FEL011046_cell_ranger_premRNA/COH073_cell_ranger/ ./FELINE_cellranger_premRNA/cell_ranger
#cp -R /home/jichen/abild/U54_data_folder/COH073/Preprocess\ FELINE\ 41-43\ premRNA\ -\ Results/qc/ FELINE_cellranger_premRNA/COH073_cell_ranger/
#cp -R /home/jichen/abild/U54_data_folder/COH073/Preprocess\ FELINE\ 41-43\ premRNA\ -\ Results/read_counts.txt.gz FELINE_cellranger_premRNA/COH073_cell_ranger/

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

