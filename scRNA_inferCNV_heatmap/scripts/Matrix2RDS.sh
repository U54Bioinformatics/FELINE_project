#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=240G
#SBATCH --time=9:00:00
#SBATCH --output=Matrix2RDS.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

export PYTHONPATH=/home/jichen/software/BETSY/install/envs/python3/lib/python3.6/site-packages/

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

FILE=`ls ./*.txt | head -n $N | tail -n 1`
prefix=${FILE%.txt}
RDS=$prefix\.RDS
echo $prefix
if [ ! -e $RDS ]; then
    echo "converting $FILE"
    /home/jichen/software/BETSY/install/envs/python3/bin/python Matrix2RDS.py --input $FILE
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

