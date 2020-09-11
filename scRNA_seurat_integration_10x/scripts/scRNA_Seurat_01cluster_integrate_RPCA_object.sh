#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --time=80:00:00
#SBATCH --output=scRNA_Seurat_01cluster_integrate_RPCA_object.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


export PATH=$PATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages:/home/jichen/software/Python_lib/lib/python2.7/site-packages/:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
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


###########################################
#input file:
#10x count file:      FEL011_10x_count01.noaggr.counts.txt
#Betsy cell metafile: FEL011_10x_cell_metadata.txt
     
#output file:
#singleR cell types:  FEL011_10x_cell_metadata.singleR_annotation.txt
#seurat project file: FEL011_10x_Seurat_2kgenes_vst_cc.rds
#seurat plots:        FEL011_10x_Seurat_2kgenes_vst_cc.pdf
###########################################

sample=FEL011046
platform=10x

if [ ! -e $sample\_$platform\_cell_metadata.UMAPcluster.txt ]; then
   echo "Running seurat to do UMAP cluster"
   #1 is to skip regression, 0 is to excute regression.
   #use 1 first to comfirm cell markers are present then use 0 to run seurat.
   #set sample to skip or keep in these two script if needed
   if [ $platform == 'icell8' ]; then
      /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_01cluster_icell8.R $sample $platform 1
   fi
   
   if [ $platform == '10x' ]; then
      echo "Standard integration"
      #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_01cluster_integrate_standard_object.R $sample $platform 0
      echo "RPCR integration"
      /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_01cluster_integrate_RPCA_object.R $sample $platform 0
      echo "SCT integration"
      #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_01cluster_integrate_SCT_object.R $sample $platform 0
      #merge only, integration of different patient or sample may not be very good due to batch effect
      #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_01cluster_merge_object.R $sample $platform 0
   fi

fi

if [ ! -e $sample\_$platform\_cell_metadata.UMAPcluster.singleR_cell_cluster_anno.txt ]; then
   echo "Running singleR to annotate cell types"
   #cell type for cells and clusters
   #FEL011_10x_cell_metadata.UMAPcluster.SingleR_anno.txt
   #FEL011_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.txt
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_02singleR_anno.R $sample $platform
fi

if [ ! -e $sample\_$platform\_seurat_figures_png ]; then
   echo "Running seurat to plot"
   #mkdir $sample\_$platform\_seurat_figures_png
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_03featureplots.R $sample $platform
   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_Seurat_03featureplots_png.R $sample $platform
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

