#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=60:00:00
#SBATCH --output=scRNA_Seurat_extract_celltype_obj.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh

#memory and time
#160G and 16 hours for singleR
#30G and 20 min for seurat

#export PATH=$PATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/
export PATH=$PATH:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
#export PYTHONPATH=$PYTHONPATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages:/home/jichen/software/Python_lib/lib/python2.7/site-packages/:/home/jichen/software/BETSY/install/lib/python2.7/site-packages
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/BETSY/install/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/lib/python2.7:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library

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

sample=FEL011046
platform=10x
celltypes=`cat celltype.list | head -n $N | tail -n 1`

if true; then
   #for celltype in Macrophages T_cells Fibroblasts Plasma_cells B_cells Adipocytes Endothelial_cells Epithelial_cells
   for celltype in $celltypes
   do
       echo "subclassify $celltype"
       if [ ! -e $sample\_$celltype\_$platform\.expression_counts.txt ]; then
           /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_extract_celltype_obj.R $sample $platform $celltype
       fi
   done
fi

############################################################################################
#Final obj:    FEL011043_T_cells_10x_Seurat_2kgenes_vst_cc.rds
#Final anno:   FEL011043_T_cells_10x.metadata.txt
############################################################################################

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

