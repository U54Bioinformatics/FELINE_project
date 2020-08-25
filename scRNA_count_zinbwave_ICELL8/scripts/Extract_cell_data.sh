#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=9:00:00
#SBATCH --output=Extract_cell_data.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./



start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=2
fi

echo "CPU: $CPU"
echo "N: $N"

#sample=FEL011043_Macrophages
#sample=FEL011043_T_cells
#sample=FEL011043_B_cells
#sample=FEL011043_Endothelial_cells
#sample=FEL011043_Adipocytes
#sample=FEL011043_Fibroblasts
#sample=FEL011043_Epithelial_cells
sample=FEL001010
platform=icell8
celltype=`cat FEL001010_icell8.celltype.txt | head -n $N | tail -n 1`
prefix=$sample\_$celltype\_$platform

#slurm
if [ ! -e $prefix\.txt ]; then
   echo "Extract meta and count data for Cell.ID provided"
   python Extract_cell_meta.py --cell $prefix\.Cell_ID.txt --meta $sample\_$platform\.metadata.txt --output $prefix\.metadata.txt
   python Extract_cell_count.py --cell $prefix\.Cell_ID.txt --count $sample\_$platform\.counts.txt --output $prefix\.count.txt
fi


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

