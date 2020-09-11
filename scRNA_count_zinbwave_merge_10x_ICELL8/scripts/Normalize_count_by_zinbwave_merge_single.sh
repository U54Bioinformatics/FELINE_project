#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=9:00:00
#SBATCH --output=Normalize_count_by_zinbwave_merge_output.sh.%A_%a.stdout
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
    N=1
fi

echo "CPU: $CPU"
echo "N: $N"

subfile=file_Epithelial_cells.txt
outfile=FEL011043_Epithelial_cells_10x
#python Normalize_count_by_zinbwave_merge_output.py --subfiles $subfile --output $outfile
python Normalize_count_by_zinbwave_merge_ssGSEA_output.py --subfiles $subfile --output $outfile


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

