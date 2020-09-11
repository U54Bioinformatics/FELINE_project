#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=FELINE_cellranger_mRNA.sh.%A_%a.stdout
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
#mkdir ./FELINE_cellranger_mRNA/COH039_FEL011_cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH039/FEL011.mRNA/FEL011_count03/ ./FELINE_cellranger_mRNA/COH039_FEL011_cell_ranger/cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH039/FEL011.mRNA/FEL011_count01.noaggr.txt ./FELINE_cellranger_mRNA/COH039_FEL011_cell_ranger/read_counts.txt

echo "COH044: FEL012,13,14,15,16"
#mkdir ./FELINE_cellranger_mRNA/COH044_FEL012_cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.mRNA/FEL012_count03/ ./FELINE_cellranger_mRNA/COH044_FEL012_cell_ranger/cell_ranger
#cp /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.mRNA/FEL012_count01.noaggr.txt ./FELINE_cellranger_mRNA/COH044_FEL012_cell_ranger/read_counts.txt
#mkdir ./FELINE_cellranger_mRNA/COH044_FEL013_cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.mRNA/FEL013_count03/ ./FELINE_cellranger_mRNA/COH044_FEL013_cell_ranger/cell_ranger
#cp /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.mRNA/FEL013_count01.noaggr.txt ./FELINE_cellranger_mRNA/COH044_FEL013_cell_ranger/read_counts.txt
#mkdir ./FELINE_cellranger_mRNA/COH044_FEL014_cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.mRNA/FEL014_count03/ ./FELINE_cellranger_mRNA/COH044_FEL014_cell_ranger/cell_ranger
#cp /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.mRNA/FEL014_count01.noaggr.txt ./FELINE_cellranger_mRNA/COH044_FEL014_cell_ranger/read_counts.txt
#mkdir ./FELINE_cellranger_mRNA/COH044_FEL015_cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.mRNA/FEL015_count03/ ./FELINE_cellranger_mRNA/COH044_FEL015_cell_ranger/cell_ranger
#cp /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.mRNA/FEL015_count01.noaggr.txt ./FELINE_cellranger_mRNA/COH044_FEL015_cell_ranger/read_counts.txt
#mkdir ./FELINE_cellranger_mRNA/COH044_FEL016_cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.mRNA/FEL016_count03/ ./FELINE_cellranger_mRNA/COH044_FEL016_cell_ranger/cell_ranger
#cp /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH044/FELINE/FEL012016.mRNA/FEL016_count01.noaggr.txt ./FELINE_cellranger_mRNA/COH044_FEL016_cell_ranger/read_counts.txt

echo "COH049: FEL017,19,20,21,22"
#mkdir ./FELINE_cellranger_mRNA/COH049_cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH049/FEL017022.mRNA/feline_p1722.noaggr.counts.txt ./FELINE_cellranger_mRNA/COH049_cell_ranger/read_counts.txt
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH049/FEL017022.mRNA/feline_p1722.10x.counts/ ./FELINE_cellranger_mRNA/COH049_cell_ranger/cell_ranger

echo "COH047: FEL023,34"
#mkdir ./FELINE_cellranger_mRNA/COH047_FEL023_cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH047/FEL023024.mRNA/FEL023_count03/ ./FELINE_cellranger_mRNA/COH047_FEL023_cell_ranger/cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH047/FEL023024.mRNA/FEL023_count01.noaggr.txt ./FELINE_cellranger_mRNA/COH047_FEL023_cell_ranger/read_counts.txt
#mkdir ./FELINE_cellranger_mRNA/COH047_FEL024_cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH047/FEL023024.mRNA/FEL024_count03/ ./FELINE_cellranger_mRNA/COH047_FEL024_cell_ranger/cell_ranger
#cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH047/FEL023024.mRNA/FEL024_count01.noaggr.txt ./FELINE_cellranger_mRNA/COH047_FEL024_cell_ranger/read_counts.txt

echo "COH055: FEL025,26,27,28"
mkdir ./FELINE_cellranger_mRNA/COH055_cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH055/FEL025028.mRNA/ ./FELINE_cellranger_mRNA/COH055_cell_ranger

echo "COH057: FEL029,30,31,32,33,34,35,36,37,38,39,40"
mkdir ./FELINE_cellranger_mRNA/COH057_cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH057/feline_p29_40_mRNA/ ./FELINE_cellranger_mRNA/COH057_cell_ranger

echo "COH064: FEL041,42,43"
mkdir ./FELINE_cellranger_mRNA/COH064_FEL041_cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH064/FEL041043.mRNA/FEL041_count03/ ./FELINE_cellranger_mRNA/COH064_FEL041_cell_ranger/cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH064/FEL041043.mRNA/FEL041_count01.noaggr.txt ./FELINE_cellranger_mRNA/COH064_FEL041_cell_ranger/read_counts.txt 
mkdir ./FELINE_cellranger_mRNA/COH064_FEL042_cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH064/FEL041043.mRNA/FEL042_count03/ ./FELINE_cellranger_mRNA/COH064_FEL042_cell_ranger/cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH064/FEL041043.mRNA/FEL042_count01.noaggr.txt ./FELINE_cellranger_mRNA/COH064_FEL042_cell_ranger/read_counts.txt
mkdir ./FELINE_cellranger_mRNA/COH064_FEL043_cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH064/FEL041043.mRNA/FEL043_count03/ ./FELINE_cellranger_mRNA/COH064_FEL043_cell_ranger/cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH064/FEL041043.mRNA/FEL043_count01.noaggr.txt ./FELINE_cellranger_mRNA/COH064_FEL043_cell_ranger/read_counts.txt

echo "COH073: FEL044,45,46"
mkdir ./FELINE_cellranger_mRNA/COH073_FEL044_cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH073/FEL044046.mRNA/FEL044_count03/ ./FELINE_cellranger_mRNA/COH073_FEL044_cell_ranger/cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH073/FEL044046.mRNA/FEL044_count03.noaggr.txt ./FELINE_cellranger_mRNA/COH073_FEL044_cell_ranger/read_counts.txt
mkdir ./FELINE_cellranger_mRNA/COH073_FEL045_cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH073/FEL044046.mRNA/FEL045_count03/ ./FELINE_cellranger_mRNA/COH073_FEL045_cell_ranger/cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH073/FEL044046.mRNA/FEL045_count03.noaggr.txt ./FELINE_cellranger_mRNA/COH073_FEL045_cell_ranger/read_counts.txt
mkdir ./FELINE_cellranger_mRNA/COH073_FEL046_cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH073/FEL044046.mRNA/FEL046_count03/ ./FELINE_cellranger_mRNA/COH073_FEL046_cell_ranger/cell_ranger
cp -R /home/jichen/Projects/Breast/scRNA/Data/Megatron_data/COH073/FEL044046.mRNA/FEL046_count03.noaggr.txt ./FELINE_cellranger_mRNA/COH073_FEL046_cell_ranger/read_counts.txt

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

