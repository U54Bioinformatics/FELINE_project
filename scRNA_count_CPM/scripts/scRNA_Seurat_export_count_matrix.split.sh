#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=350G
#SBATCH --time=90:00:00
#SBATCH --output=scRNA_Seurat_export_count_matrix.split.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


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
    N=2
fi

echo "CPU: $CPU"
echo "N: $N"

sample=FEL011046
platform=10x
patient=`cat FEL011046_10x.patients.list | head -n $N | tail -n 1`
prefix=$sample\_$platform

#extract CPM for each patient
if [ ! -e $sample\_$platform\_$patient\_gene_symbols.raw.counts.txt ]; then
   echo "Export normalized,scaled count, and CPM data"
   /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_export_count_matrix.split.R $sample $platform $patient
fi

#merge all patient into single file
if [ ! -e $sample\_$platform\_gene_symbols.raw.counts.txt ]; then
   echo "merge count"
   #python scRNA_Seurat_export_count_matrix.split.merge.py --list $prefix\.patients.list --prefix $prefix
fi


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

