#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=280G
#SBATCH --time=9:00:00
#SBATCH --output=Normalize_count_by_zinbwave_merge_multiple.sh.%A_%a.stdout
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

sample=FEL001046
platform=scRNA
celltype=`cat FEL001046.celltype.txt | head -n $N | tail -n 1`

#subfile=file_Epithelial_cells.txt
#outfile=FEL011043_Epithelial_cells_10x
subfile=file_$celltype\.txt
outfile=$sample\_$celltype\_$platform

#data='zinbwave'
data='ssGSEA'
if [ $data == 'zinbwave' ] && true; then
    if [ ! -e $outfile\.zinbwave.normalized.txt ]; then
        echo "merge zinbwave count"
        python Normalize_count_by_zinbwave_merge_zinbwave_output.py --subfiles $subfile --output $outfile
    fi
fi

if [ $data == 'ssGSEA' ] && true; then
    if [ ! -e $outfile\.zinbwave.normalized.ssGSEA_scores.txt ]; then
        echo "merge ssGSEA"
        python Normalize_count_by_zinbwave_merge_ssGSEA_output.py --subfiles $subfile --output $outfile
    fi
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

