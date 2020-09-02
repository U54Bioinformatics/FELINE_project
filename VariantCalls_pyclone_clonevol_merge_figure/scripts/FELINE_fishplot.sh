#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --output=FELINE_fishplot.sh.%A_%a.stdout
#SBATCH -p abild,all
#SBATCH --workdir=./


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

export R_LIBS=/home/jichen/software/BETSY/install/envs/R_combine_pdf_figures/lib/R/library/

#echo "sort patients"
#/home/jichen/software/BETSY/install/envs/R_combine_pdf_figures/bin/Rscript Run_merge_tree_helper.R

echo "prepare script"

# sort by respone then by arm
prefix=FELINE_fishplot_by_arm_response
if [ ! -e $prefix\.png ] && true; then
python FELINE_Step5_merge_fishplot.py --input FELINE_patient_info.sorted.revised.txt --pyclone_cnv FELINE_pyclone_model.txt --project $prefix
/home/jichen/software/BETSY/install/envs/R_combine_pdf_figures/bin/Rscript $prefix\.R
convert -density 100 -quality 100 $prefix\.pdf $prefix\.png
fi

# sort by arm then by shannon delta
# FELINE_patient_info.sorted_by_Arm_shannon_delta.revised.txt
for analysis in Arm_deltaDominance Arm_Dominance0 Arm_shannon_delta Arm_shannon_H0 deltaDominance Dominance0 shannon_delta shannon_H0
#for analysis in Arm_shannon_delta Arm_shannon_H0
do
   #analysis=Arm_shannon_delta
   prefix=FELINE_fishplot_by_$analysis
   info=FELINE_patient_info.sorted_by_$analysis\.revised.txt
   echo $analysis
   echo $prefix
   echo $info
   if [ ! -e $prefix\.png ] && true; then
       echo " Do fishplot ......"
       python FELINE_Step5_merge_fishplot.py --input $info --pyclone_cnv FELINE_pyclone_model.txt --project $prefix
       /home/jichen/software/BETSY/install/envs/R_combine_pdf_figures/bin/Rscript $prefix\.R
       convert -density 100 -quality 100 $prefix\.pdf $prefix\.png
   fi
done



end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

