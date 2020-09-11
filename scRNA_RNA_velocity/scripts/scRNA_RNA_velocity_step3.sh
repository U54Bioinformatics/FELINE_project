#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=280G
#SBATCH --time=40:00:00
#SBATCH --output=scRNA_RNA_velocity_step3.sh.%A_%a.stdout
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

sample=`cat FELINE_sample.list | head -n $N | tail -n 1`
#sample=FEL011046_ArmA
#sample=FEL011046_ArmB
#sample=FEL011046_ArmC
#sample=FELINE_test_ArmA
#sample=FELINE_test_ArmB
loom_file=$sample\.loom
# spliced or unspliced
spliced=unspliced
if [ $sample == 'FEL011046_ArmA' ]; then
   ln -s FEL011046.loom $loom_file
   patients="FEL014 FEL015 FEL016 FEL018 FEL020 FEL023 FEL028 FEL029 FEL041 FEL042 FEL043 FEL044 FEL045"
   umap_pos=Jason_input/10_umap_visualization/ARM_A_cancer_cell_UMAP_results.umap.txt
   gene_list=Jason_input/10_umap_visualization/ARM_A_cancer_cell_UMAP_results.gene.txt
fi
if [ $sample == 'FEL011046_ArmB' ]; then
   ln -s FEL011046.loom $loom_file
   patients="FEL011 FEL012 FEL017 FEL019 FEL022 FEL025 FEL027 FEL032 FEL033 FEL034 FEL037 FEL046"
   umap_pos=Jason_input/10_umap_visualization/ARM_B_Day180_cancer_cell_UMAP_results.umap.txt
   gene_list=Jason_input/10_umap_visualization/ARM_B_Day180_cancer_cell_UMAP_results.gene.txt
fi
if [ $sample == 'FEL011046_ArmC' ]; then
   ln -s FEL011046.loom $loom_file
   patients="FEL013 FEL021 FEL024 FEL026 FEL030 FEL031 FEL035 FEL036 FEL038 FEL039 FEL040"
   umap_pos=Jason_input/10_umap_visualization/ARM_C_Day180_cancer_cell_UMAP_results.umap.txt
   gene_list=Jason_input/10_umap_visualization/ARM_C_Day180_cancer_cell_UMAP_results.gene.txt
fi
if [ $sample == 'FELINE_test_ArmA' ]; then
   patients="FEL014 FEL015 FEL016"
   ln -s FELINE_test.loom $loom_file
   umap_pos=Jason_input/10_umap_visualization/ARM_A_cancer_cell_UMAP_results.umap.txt
   gene_list=Jason_input/10_umap_visualization/ARM_A_cancer_cell_UMAP_results.gene.txt
fi
if [ $sample == 'FELINE_test_ArmB' ]; then
   ln -s FELINE_test.loom $loom_file
   patients="FEL011 FEL012 FEL017"
   umap_pos=Jason_input/10_umap_visualization/ARM_B_Day180_cancer_cell_UMAP_results.umap.txt
   gene_list=Jason_input/10_umap_visualization/ARM_B_Day180_cancer_cell_UMAP_results.gene.txt
fi

# this seurat obj can be a pre computed obj with UMAP coordinate (set umap_txt=0)
# I use a txt umap coordinate from Jason (set umap_txt=1)
obj_velocity=$sample\_velocity.rds

echo "Step3: Analyze RNA velocity and project the results on UMAP"
echo "$sample, $umap_pos, $gene_list, $spliced"
#export PATH=$PATH:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
#export PYTHONPATH=$PYTHONPATH:/home/jichen/software/BETSY/install/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/lib/python2.7:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_velocityR/lib/R/library
module load hdf5/1.10.3

denovo_umap=0
umap_txt=1
#if [ $denovo_umap == 1 ]; then
#   echo "Perform RNA velocity on denovo umap"
#   #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_RNA_velocity_run_velocityR_denovo_umap.R $sample $loom_file
#   /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_RNA_velocity_run_velocityR_denovo_umap_subset_genes.R $sample $loom_file
#fi
if [ $denovo_umap == 0 ]; then
   if [ $umap_txt == 0 ]; then
       echo "Perform RNA velocity on pre estimated umap in seurat obj"
       /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_RNA_velocity_run_velocityR_given_umap.R $sample $loom_file $obj_seurat
   else
       echo "Perform RNA velocity on pre estimated umap in txt file"
       # plot umap in txt
       #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_RNA_velocity_plot_umap_in_txt.R $sample $loom_file $umap_pos
       # all cells and all genes
       /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_RNA_velocity_run_velocityR_given_umap_in_txt.R $sample $loom_file $umap_pos
       # all cells and only genes used to create umap (failed)
       #/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_RNA_velocity_run_velocityR_given_umap_in_txt_subset_genes.R $sample $loom_file $umap_pos $gene_list $spliced
       # cells from one patient and all genes
       #echo $patients
       #for patient in $patients
       #do
       #    echo "processing $patient"
       #    /home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript scRNA_RNA_velocity_run_velocityR_given_umap_in_txt_patient.R $sample $loom_file $umap_pos $patient
       #done
   fi
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

