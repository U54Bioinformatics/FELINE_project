#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --time=9:00:00
#SBATCH --output=infercnv_heatmap_individual_subclone.sh.%A_%a.stdout
#SBATCH -p fast,all,abild
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
    N=3
fi

echo "CPU: $CPU"
echo "N: $N"

platform=10x
#patient=FEL028P131
patient=FEL045P122
#patient=FEL026P140
#patient=FEL025P142
#patient=`cat patients.list | head -n $N | tail -n 1`

if [ ! -e $patient\_HMM_infersubclone_k$N\.png ]; then
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
   ~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript infercnv_heatmap_individual_plot_sample_all_subclone.R $patient 0 $N $infer_clust
   #https://github.com/kundajelab/phantompeakqualtools/issues/3
   #~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript --max-ppsize=500000 infercnv_heatmap_individual_plot_sample_all.R $patient 0 2 $infer_clust
   #cat infercnv_heatmap_individual_plot_sample_all.R | R --slave
fi
   
 
#plot infercnv final obversation based on HMM clustering rank
#infer_noclust="final"
#~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript infercnv_heatmap_individual_plot_sample_rank_by_HMM.R $patient 1 0 $infer_noclust
infer_noclust="final"
~/software/Miniconda2/miniconda2/envs/Infercnv_heatmap/bin/Rscript infercnv_heatmap_individual_plot_sample_all_subclone_rank_by_HMM.R $patient 1 $N $infer_noclust

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

