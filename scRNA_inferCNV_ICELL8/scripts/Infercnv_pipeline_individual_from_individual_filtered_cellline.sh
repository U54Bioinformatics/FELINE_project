#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=220G
#SBATCH --time=90:00:00
#SBATCH --output=Infercnv_pipeline_individual_from_individual_filtered_cellline.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


export PATH=$PATH:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/BETSY/install/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/lib/python2.7:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
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


##############################################
#input
#metadata, two cols with Cell and Sample
#Cell    Sample
#FEL011_FEL011_M_AAGTACCGTCCGGCAT        FEL011_M
#FEL011_FEL011_S_ACGGTCGGTTGCTCAA        FEL011_S

#############################################



#sample=FEL001016
platform=10x
downsample=all
# cellline: if control is "cellline" use cell line/other patients normal as control, need to merge cell line data into the matrix
# normal  : otherwise use normal from patients themselve as control (need to run after singleR annoation), just copy the data
control="cellline" #cellline, normal
#control_pre="FELINE.Normal"
control_pre="FELINE.icell8.Immune" #HMEC, Fibroblasts, 97OVCZ.ImmuneCells (Immune), GTEx (Breast_normal) Fibroblasts_HMEC_97OVCZ (edit in infercnv_plot.R)
control_cell="Immune"
control_count=~/Projects/Breast/scRNA/Data/scRNA_cell_lines/$control_pre\.counts.txt
control_annot=~/Projects/Breast/scRNA/Data/scRNA_cell_lines/$control_pre\.annot.txt
patient=`cat patients.list | head -n $N | tail -n 1`
###########################################################

###########################################################

#create raw matrix and raw annotation files from 10x datadir
if [ ! -e $patient\_$platform\_counts.raw.matrix ] && [ $control == "cellline" ]; then
   echo "Create raw matrix"
   ##make sure metadata is "Cell\tSample" and remove normal cell if needed
   echo "Cell"$'\t'"Sample" > $patient\_$platform\_cell_metadata.txt
   #do only epithelial cell
   #awk -F"\t" '$25~/Epithelial/' $patient\_$platform\_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt | grep -v "^Cell" | cut -f1,5 | shuf >> $patient\_$platform\_cell_metadata.txt
   #do all cell types
   cat $patient\_$platform\_cell_metadata.UMAPcluster.txt | grep -v "^Cell" | cut -f1,5 | shuf >> $patient\_$platform\_cell_metadata.txt
   if [ $platform == '10x' ]; then
      #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript Infercnv_convert_to_infercnv_input_10x.R $patient $platform $control_count $control_annot
      /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript Infercnv_convert_to_infercnv_input_10x_filtered.R $patient $platform $control_count $control_annot
   fi
   if [ $platform == 'icell8' ]; then
      /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript Infercnv_convert_to_infercnv_input_icell8.R $patient $platform $control_count $control_annot
   fi
fi

if [ ! -e $patient\_$platform\_counts.raw.matrix ] && [ $control == "normal" ]; then
   echo "Copying raw.matrix file from filtered.counts.txt"
   cp $patient\_$platform\_gene_symbols.filtered.counts.txt $patient\_$platform\_counts.raw.matrix
   cut -f1,24 $patient\_$platform\_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt | grep -v "^Cell" | shuf > $patient\_$platform\_counts.anno.raw.txt
fi

#match genes from 10x and hg19 gene position file
if [ ! -e $patient\_$platform\_counts.matrix ]; then
   echo "Make gene position file"
   python Infercnv_match_matrix_and_gene_pos_files.py --matrix $patient\_$platform\_counts.raw.matrix
fi

#format anno file and downsample if needed
if [ ! -e $patient\_$platform\_counts.anno.$downsample\.txt ]; then
   echo "Make annotation file"
   python Infercnv_convert_anno_input_down.py --input $patient\_$platform\_counts.anno.raw.txt --downsample $downsample > $patient\_$platform\_counts.anno.$downsample\.txt
fi

#export PKG_CONFIG_PATH=/home/jichen/software/JAGS/JAGS/lib/pkgconfig/
#export LD_RUN_PATH=/home/jichen/software/JAGS/JAGS/lib/JAGS/
#export R_LIBS=/home/jichen/software/R.3.6.0/R_package_3.6.0:$R_LIBS
export PKG_CONFIG_PATH=/home/jichen/software/Miniconda2/miniconda2/envs/Infercnv/lib/pkgconfig
export LD_RUN_PATH=/home/jichen/software/Miniconda2/miniconda2/envs/Infercnv/lib/JAGS/
export R_LIBS=$R_LIBS:/home/jichen/software/Miniconda2/miniconda2/envs/Infercnv/lib/R/library/
if [ true == true ]; then
   echo "Run infercnv: about ~6h on apollo"
   #ln -s $sample\_$platform\_counts.matrix $patient\_$platform\_counts.matrix
   #ln -s $sample\_$platform\_counts.gene_pos.txt $patient\_$platform\_counts.gene_pos.txt
   #if [ ! $downsample == 'all' ]; then 
   #   echo "Run infercnv for subset cells"
   #   /home/jichen/software/Miniconda2/miniconda2/envs/Infercnv/bin/Rscript Infercnv_pipeline_individual_from_individual_filtered_cellline_plot.R $patient $platform $downsample $control_cell
   #   python Infercnv_rename_files.py --infercnv_outdir $patient\_$platform\_counts.$downsample\.output_dir_plot --prefix $patient\_$platform\_$downsample\_$control
   #fi
   #if [ $downsample == 'all' ]; then
   
      echo "Run infercnv for all cells (no plot from infercnv, faster)"
      #min_cells_per_gene 
      cells=100
      /home/jichen/software/Miniconda2/miniconda2/envs/Infercnv/bin/Rscript Infercnv_pipeline_individual_from_individual_filtered_cellline_plot_all.R $patient $platform $downsample $control_cell $cells
      #python Infercnv_rename_files.py --infercnv_outdir $patient\_$platform\_counts.$downsample\.output_dir_plot_$cells"_cells" --prefix $patient\_$platform\_$downsample\_$control
      #/home/jichen/software/Miniconda2/miniconda2/envs/Infercnv/bin/Rscript Infercnv_CNA_correlation.R $patient final
      #/home/jichen/software/Miniconda2/miniconda2/envs/Infercnv/bin/Rscript Infercnv_CNA_correlation.R $patient hmm
   #fi
fi


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

