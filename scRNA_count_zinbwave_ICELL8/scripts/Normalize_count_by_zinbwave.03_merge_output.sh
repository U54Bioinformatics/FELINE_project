#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=280G
#SBATCH --time=9:00:00
#SBATCH --output=Normalize_count_by_zinbwave.03_merge_output.sh.%A_%a.stdout
#SBATCH -p all,abild
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

sample=FEL001010
platform=icell8
celltype=`cat FEL001010_icell8.celltype.txt | head -n $N | tail -n 1`
subfile=file_$celltype\.txt
outfile=$sample\_$celltype\_$platform
#subfile=file_T_cells.txt
#outfile=FEL011043_T_cells_10x
#subfile=file_Fibroblasts.txt
#outfile=FEL011043_Fibroblasts_10x
#subfile=file_Epithelial_cells.txt
#outfile=FEL011043_Epithelial_cells_10x
#subfile=file_Endothelial_cells.txt
#outfile=FEL011043_Endothelial_cells_10x
python Normalize_count_by_zinbwave_merge_output.py --subfiles $subfile --output $outfile
#python Normalize_count_by_zinbwave_merge_ssGSEA_output.py --subfiles $subfile --output $outfile


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

