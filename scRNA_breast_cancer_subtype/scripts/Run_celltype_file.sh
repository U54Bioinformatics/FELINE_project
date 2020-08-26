#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --output=Run_celltype_file.sh.%A_%a.stdout
#SBATCH -p fast,all,abild
#SBATCH --workdir=./


start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=3
fi

echo "CPU: $CPU"
echo "N: $N"

export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library

/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript Run_celltype_file.R
/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript Run_gene_id_file.R


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

