#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --output=Run_genefu_scRNA.sh.%A_%a.stdout
#SBATCH -p fast,all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh



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

export R_LIBS=/home/jichen/software/BETSY/install/envs/breast_cancer_subtype/lib/R/library

cpm_file=`ls *.CPM.txt | head -n $N | tail -n 1`
prefix=${cpm_file%.CPM.txt}
celltype=FEL011046.cell_type.txt
geneid=FEL011046.gene_id.txt
echo $prefix

if [ ! -e $prefix\.genefu.sum.txt ]; then
   echo "Run genefu: $prefix"
   /home/jichen/software/BETSY/install/envs/breast_cancer_subtype/bin/Rscript Run_genefu_scRNA.R $cpm_file $celltype $geneid
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

