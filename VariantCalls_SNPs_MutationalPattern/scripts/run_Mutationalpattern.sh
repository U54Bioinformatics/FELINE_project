#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --output=run_Mutationalpattern.sh.%A_%a.stdout
#SBATCH -p all
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


#export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_mutational_pattern/lib/R/library:$R_LIBS
export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_mutational_pattern_R3.6/lib/R/library

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


/home/jichen/software/BETSY/install/envs/DNA_mutational_pattern/bin/Rscript Mutationalpattern_FELINE.R


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

