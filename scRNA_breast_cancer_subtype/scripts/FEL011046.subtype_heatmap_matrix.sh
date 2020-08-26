#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=4:00:00
#SBATCH --output=FEL011046.subtype_heatmap_matrix.sh.%A_%a.stdout
#SBATCH -p fast,all
#SBATCH --workdir=./

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

#############################################################################################################
#############################################################################################################

prefix=FEL011046.subtype_heatmap_matrix
echo "prefix: $prefix"

#export R_LIBS=/home/jichen/software/Miniconda2/miniconda2/envs/ggplot/lib/R/library/
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library

if [ ! -e $prefix\.heatmap\.png ]; then
    echo "plot dotplot and heatmap"
    /home/jichen/software/Miniconda2/miniconda2/envs/ggplot/bin/Rscript FEL011046.subtype_heatmap_matrix.R $prefix
    convert -density 300 -quality 100 $prefix\.heatmap\.pdf $prefix\.heatmap\.png
fi

echo "Done"

