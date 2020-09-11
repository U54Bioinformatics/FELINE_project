#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=9:00:00
#SBATCH --output=FEL001010_scRNA.cellcycle_barplot.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library

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

/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript FEL001010_scRNA.celltype_barplot.R


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

