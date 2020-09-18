#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=9:00:00
#SBATCH --output=scRNA_Seurat_count2umap.sh.%A_%a.stdout
#SBATCH -p abild,all
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


export PATH=$PATH:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/BETSY/install/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/lib/python2.7:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
#export PYTHONPATH=/home/jichen/software/BETSY/install/envs/scRNA_seurat_wrapper/lib/python3.7/site-packages/
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
#export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_seurat_wrapper/lib/R/library

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

###########################################
#input file:
#10x count file:      FEL011_10x_count01.noaggr.counts.txt
#Betsy cell metafile: FEL011_10x_cell_metadata.txt
     
#output file:
#seurat project file: FEL011_10x_Seurat_2kgenes_vst_cc.rds
###########################################

#sample=FEL012016
#sample=FEL023P121
#platform=10x
#sample=`cat patients.list | head -n $N | tail -n 1`
#prefix=FEL011040_chunk0_10x_gene_symbols.filtered.counts.T_cells
#prefix=FEL011040_10x_gene_symbols.filtered.counts.Epithelial_cells
#prefix=FEL011043_Epithelial_cells_10x
#prefix=FEL011043_Adipocytes_10x
prefix=FEL001010_icell8.immune
method=Seurat
#method=FastMNN
R=/home/jichen/software/BETSY/install/envs/scRNA_seurat_wrapper/bin/Rscript
if [ ! -e $prefix\_filtered.counts_png ]; then
    mkdir $prefix\_filtered.counts_png
    $R scRNA_Seurat_count2umap.R $prefix filtered.counts $method 1
fi


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

