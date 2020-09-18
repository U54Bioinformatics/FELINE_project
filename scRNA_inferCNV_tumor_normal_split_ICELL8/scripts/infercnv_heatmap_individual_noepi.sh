#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=80:00:00
#SBATCH --output=infercnv_heatmap_individual_noepi.sh.%A_%a.stdout
#SBATCH -p abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh

export R_LIBS=/home/jichen/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/lib/R/library/
source activate Infercnv_heatmap
#fix Error: node stack overflow
ulimit -s 500000
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

platform=10x
patient=`cat patients.list | head -n $N | tail -n 1`

if [ ! -e $patient\_CNA_infer_NoClust\.png ]; then
   echo "Plotting $patient ..."
   #head -n 1 cell_metadata.txt > $patient\_cell_metadata.txt
   #grep $patient cell_metadata.txt >> $patient\_cell_metadata.txt
   #parameters
   #1. patient id: FEL011
   #2. nocluster:  0 to skip nocluster, otherwise plot nocluster
   #3. cluster number: 0 to skil cluster, otherwise do clustering using the number
   #4. infercnv observation file type: pre (preliminary), hmm (hmm14), final
   infer_noclust="final"
   infer_clust="hmm"
   echo "plot 2 clusters using hmm to split tumor and normal"
   #~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all.R $patient 0 2 $infer_clust
   #~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all.R $patient 0 3 $infer_clust
   #~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all.R $patient 0 4 $infer_clust
   #~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all.R $patient 0 5 $infer_clust
   ~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all.R $patient 0 6 $infer_clust
   #~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all.R $patient 0 8 $infer_clust
   echo "plot 2 clusters using final to double check hmm clusters"
   #~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all_rank_by_HMM.R $patient 1 2 $infer_noclust
   #~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all_rank_by_HMM.R $patient 1 3 $infer_noclust
   #~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all_rank_by_HMM.R $patient 1 4 $infer_noclust 
   #~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all_rank_by_HMM.R $patient 1 5 $infer_noclust
   ~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all_rank_by_HMM.R $patient 1 6 $infer_noclust
   #~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all_rank_by_HMM.R $patient 1 8 $infer_noclust
   #cat infercnv_heatmap_individual_plot_sample_all.R | R --slave
fi
   


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

