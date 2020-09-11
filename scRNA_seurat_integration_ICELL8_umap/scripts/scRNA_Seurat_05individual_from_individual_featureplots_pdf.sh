#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=280G
#SBATCH --time=9:00:00
#SBATCH --output=scRNA_Seurat_05individual_from_individual_featureplots_pdf.sh.%A_%a.stdout
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

R=/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript

#sample=FEL011043_Fibroblasts
#sample=FEL011043_Epithelial_cells
sample=FEL001010_integrated_RPCA
#sample=FEL011043_integrated_standard
#sample=FEL011043_integrated_SCT
#sample=FEL011043_integrated_RPCA_default
platform=icell8

if [ ! -e $sample\_$platform\_Seurat_UMAP_Celltype ]; then
    mkdir $sample\_$platform\_Seurat_UMAP_Celltype
    # all cell type
    #$R scRNA_Seurat_05individual_from_individual_featureplots_pdf.R $sample $platform
    # four broad cell type: cancer/normal/immune/stromal 
    $R scRNA_Seurat_05individual_from_individual_featureplots_pdf_broad_class.R $sample $platform

fi


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

