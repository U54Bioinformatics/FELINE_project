#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=250G
#SBATCH --time=40:00:00
#SBATCH --output=Run_MEDALT.sh.%A_%a.stdout
#SBATCH -p abild,all
#SBATCH --workdir=./


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

export R_LIBS=/home/jichen/software/BETSY/install/envs/MEDALT_sc_tree/lib/R/library
export PATH=/home/jichen/software/BETSY/install/envs/MEDALT_sc_tree/bin/:$PATH
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/MEDALT_sc_tree/lib/python2.7/site-packages/
MEDALT=/home/jichen/software/Single_cell_RNAseq/MEDALT/MEDALT-master/

#cnvfile=./example/FEL015_500ran.CNV.txt
#outdir=./example/outputRNAT_500
#cnvfile=./example/scRNA.CNV.txt
#outdir=./example/outputRNAT
#python scTree.py -P ./ -I $cnvfile -O $outdir -D R -G hg38

patient=`cat patients.list | head -n $N | tail -n 1`

Infercnv_obs_origin_dir=FEL011046_infercnv_HMMi6.observations
Infercnv_obs_MEDALT_dir=FEL011046_infercnv_HMMi6.observations_4MEDALT
Infercnv_subclone_dir=FEL011046_infercnv_subclones_final
Infercnv_MEDALT_result_dir=FEL011046_infercnv_MEDALT_results
gene_expr_rds=FEL011046_data_gene_pathway.RDS



cnv_file_raw=$Infercnv_obs_origin_dir/$patient\.infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt
cnv_file_new=$Infercnv_obs_MEDALT_dir/$patient\.HMMi6.CNV.txt
mkdir $Infercnv_obs_MEDALT_dir
if [ ! -e $cnv_file_new ]; then
   echo "Step1: Prepare infercnv data for MEDALT by removing quote and replace space with tab"
   sed s'/"//g' $cnv_file_raw | sed s'/\s/\t/g' > $cnv_file_new
fi

cnv_file_rand=$Infercnv_MEDALT_result_dir/$patient\.HMMi6.rand500.CNV.txt
if [ ! -e $cnv_file_rand ]; then
   echo "Step2: Prepare infercnv data for final cells and final rand500 cells"
   /home/jichen/software/BETSY/install/envs/MEDALT_sc_tree/bin/Rscript Run_MEDALT_prepare_cell_list.R $patient $Infercnv_subclone_dir $Infercnv_obs_MEDALT_dir $Infercnv_MEDALT_result_dir $gene_expr_rds
fi


# rand500: relatively faster
cnvfile=./$Infercnv_MEDALT_result_dir/$patient\.HMMi6.rand500.CNV.txt
## for patients with fewer than 500 cells, link the files.
if [ ! -e $cnvfile ]; then
   mkdir $Infercnv_MEDALT_result_dir
   echo "Step3: Use final files if patients have fewer than 500 cells"
   ln -s `pwd`/$Infercnv_MEDALT_result_dir/$patient\.HMMi6.final.CNV.txt `pwd`/$Infercnv_MEDALT_result_dir/$patient\.HMMi6.rand500.CNV.txt
   ln -s `pwd`/$Infercnv_MEDALT_result_dir/$patient\.HMMi6.final.CNV.names.txt `pwd`/$Infercnv_MEDALT_result_dir/$patient\.HMMi6.rand500.CNV.names.txt
fi
outdir=./$Infercnv_MEDALT_result_dir/$patient\_outputRNAT_rand500
#cnvfile=./FEL015_500ran.CNV.txt
#outdir=./outputRNAT_rand500
if [ ! -e $outdir/$patient\.CNV.tree.txt ] && false; then
   echo "Step3: Analysis scRNA tree with MEDALT for rand500 final cells"
   echo "$outdir"
   mkdir $outdir
   python $MEDALT/scTree.py -P $MEDALT -I $cnvfile -O $outdir -D R -G hg38
   mv $outdir/CNV.tree.txt $outdir/$patient\.CNV.tree.txt
fi

# all cells: can be very slow and have memory issues
cnvfile_all=./$Infercnv_MEDALT_result_dir/$patient\.HMMi6.final.CNV.txt
outdir_all=./$Infercnv_MEDALT_result_dir/$patient\_outputRNAT_final
if [ ! -e $outdir_all/$patient\.CNV.tree.txt ] && false; then
   echo "Step3: Analysis scRNA tree with MEDALT for all final cells"
   echo "$outdir_all"
   mkdir $outdir_all
   python $MEDALT/scTree.py -P $MEDALT -I $cnvfile_all -O $outdir_all -D R -G hg38
   mv $outdir_all/CNV.tree.txt $outdir_all/$patient\.CNV.tree.txt
fi

# plot scRNA MEDALT tree with igraph using 500 cells
# change false to true to run
figure=$outdir/$patient\.singlecell.tree.sample_color.pdf
if [ ! -e $figure ] && false; then
   echo "Step4: Plot scRNA with colors annotated"
   /home/jichen/software/BETSY/install/envs/MEDALT_sc_tree/bin/Rscript Run_MEDALT_plot_sc_tree.R $patient $Infercnv_subclone_dir $outdir $gene_expr_rds
fi

# plot scRNA MEDALT tree with igraph using all cells
# change false to true to run
figure_all=$outdir_all/$patient\.singlecell.tree.sample_color.pdf
if [ ! -e $figure_all ] && false; then
   echo "Step4: Plot scRNA with colors annotated for all final cells"
   /home/jichen/software/BETSY/install/envs/MEDALT_sc_tree/bin/Rscript Run_MEDALT_plot_sc_tree.R $patient $Infercnv_subclone_dir $outdir_all $gene_expr_rds
fi


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

