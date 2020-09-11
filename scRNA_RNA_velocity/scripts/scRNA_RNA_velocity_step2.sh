#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=9:00:00
#SBATCH --output=scRNA_RNA_velocity_step2.sh.%A_%a.stdout
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

export PYTHONPATH=/home/jichen/software/BETSY/install/envs/scRNA_velocyto/lib/python3.6/site-packages/
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_velocityR/lib/R/library
module load hdf5/1.10.3

# Run FEL011046
project=FEL011046
sample=samples.list
# Run subset test
#project=FELINE_test
#sample=samples.test.list
bam_dir=FEL011046_cell_ranger.mRNA

echo "Step2: If there are multiple loom from different samples/patients do some merging here to generate a single loom that is cresponding to the UMAP to analyze"
#The R loomR failed to merge multiple loom files
#/home/jichen/software/BETSY/install/envs/scRNA_velocityR/bin/Rscript scRNA_RNA_velocity_merge_loom_files.R $bam_dir $sample $project\.loom
/home/jichen/software/BETSY/install/envs/scRNA_velocyto/bin/python scRNA_RNA_velocity_merge_loom_files.py --bamdir $bam_dir --sample $sample --output $project

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

