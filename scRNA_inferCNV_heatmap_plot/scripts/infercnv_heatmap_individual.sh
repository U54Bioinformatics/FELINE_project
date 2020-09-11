#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=9:00:00
#SBATCH --output=infercnv_heatmap_individual.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh

export R_LIBS=/home/jichen/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/lib/R/library/
source activate Infercnv_heatmap
ulimit -s 36384
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


   #parameters
   #1. patient id: FEL011
   #2. nocluster:  0 to skip nocluster, otherwise plot nocluster
   #3. cluster number: 0 to skil cluster, otherwise do clustering using the number
   #4. infercnv observation file type: pre (preliminary), hmm (hmm14), final
   infer_noclust="final"
   infer_clust="hmm"
   platform=10x
   patient=FEL001010
   if [ ! -e $patient\_CNA_infer_NoClust\.png ]; then
       ~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript infercnv_heatmap_individual_plot_sample.R $patient 1 0 $infer_noclust
   fi
   patient=FEL011046
   if [ ! -e $patient\_CNA_infer_NoClust\.png ]; then
       ~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript infercnv_heatmap_individual_plot_sample.R $patient 1 0 $infer_noclust
   fi
   

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

