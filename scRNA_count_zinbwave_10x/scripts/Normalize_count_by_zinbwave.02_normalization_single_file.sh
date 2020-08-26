#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=50G
#SBATCH --time=9:00:00
#SBATCH --output=Normalize_count_by_zinbwave.02_normalization_single_file.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

export R_LIBS=/home/jichen/software/BETSY/install/envs/zinbwave/lib/R/library
export OMP_NUM_THREADS=1
export USE_SIMPLE_THREADED_LEVEL3= 1

start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=6
fi

echo "CPU: $CPU"
echo "N: $N"

#count=FEL011040_10x_gene_symbols.filtered.counts.Adipocytes.txt
#count=FEL011040_10x_gene_symbols.filtered.counts.B_cells.txt
#count=FEL011040_10x_gene_symbols.filtered.counts.Unclassified.txt
#count=`cat file.list | head -n $N | tail -n 1`
#/home/jichen/software/BETSY/install/envs/zinbwave/bin/Rscript Normalize_count_by_zinbwave.R $count 


celltype=`cat FEL011046_10x.celltype.txt | head -n $N | tail -n 1`
echo "Split files"
count=`cat file_$celltype\.txt | head -n 1 | tail -n 1`
#count=`cat file_T_cells.txt | head -n $N | tail -n 1`
#count=`cat file_Macrophages.txt | head -n $N | tail -n 1`
#count=`cat file_Epithelial_cells.txt | head -n $N | tail -n 1`
#count=`cat file_Fibroblasts.txt | head -n $N | tail -n 1`
#count=`cat file_Adipocytes.txt | head -n $N | tail -n 1`
#count=`cat file_Endothelial_cells.txt | head -n $N | tail -n 1`

prefix=${count%.txt}
echo $count
echo $prefix
if [ ! -e $prefix\.zinbwave.normalized.txt ]; then
   echo 'Run zinbwave'
   /home/jichen/software/BETSY/install/envs/zinbwave/bin/Rscript Normalize_count_by_zinbwave.split.R $count
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

