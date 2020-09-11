#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=8:00:00
#SBATCH --output=cp_files.sh.%A_%a.stdout
#SBATCH -p fast,all,abild
#SBATCH --workdir=./

#sbatch --array 1-3 cp_files.sh

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

#FILE=`cat file_bai.list | head -n $N | tail -n 1`
FILE=`cat file_bam.list | head -n $N | tail -n 1`
NAME=`basename $FILE`

echo $FILE
echo $NAME

if [ ! -e $NAME ]; then
echo "copy $FILE ......"
 
   cp $FILE ./
   #rsync -v -L $FILE ./
   #rsync -v -e "sshpass -f '/home/jichen/software/rsync/tikavah.txt' ssh -o StrictHostKeyChecking=no -p 22" --progress jichen@10.30.83.100:$FILE ./
 
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

