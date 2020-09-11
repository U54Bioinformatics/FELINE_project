#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --time=9:00:00
#SBATCH --output=FEL011046_data_gene_pathway.sh.%A_%a.stdout
#SBATCH -p abild,all
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

cat FEL011046_data_gene_pathway.R | /home/jichen/software/BETSY/install/envs/scRNA/bin/R --slave

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

