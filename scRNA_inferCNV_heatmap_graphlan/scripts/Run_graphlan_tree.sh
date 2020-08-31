#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --output=Run_graphlan_tree.sh.%A_%a.stdout
#SBATCH -p abild,all
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

export PATH=/home/jichen/software/BETSY/install/envs/graphlan/bin/:$PATH
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/graphlan/lib/python2.7/site-packages/

patient=`cat patients_arm_response.list | cut -f1 | head -n $N | tail -n 1`
#annot=annot_0.txt
annot=$patient\.annot1.txt

graphlan_annotate.py --annot $annot MEDALT_tree/$patient\.tree $patient\.tree.xml

#pdf/svg are not vector, ps works as verctor but hard to edit
#graphlan.py $patient\.tree.xml $patient\.tree.ps --format ps --size 3.5 --pad 0.05 --avoid_reordering
#png
graphlan.py $patient\.tree.xml $patient\.tree.png --dpi 300 --size 3.5 --pad 0.05 --avoid_reordering


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

