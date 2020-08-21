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
prefix=FELINE_fishplot
# Non-responder
if true; then
python FELINE_Step5_merge_fishplot.py --input FELINE_patient_info.sorted.revised.txt --pyclone_cnv FELINE_pyclone_model.txt --project $prefix
/home/jichen/software/BETSY/install/envs/R_combine_pdf_figures/bin/Rscript $prefix\.R
convert -density 100 -quality 100 $prefix\.pdf $prefix\.png
#mv $prefix\.pdf $prefix\.Non-responder.pdf
#mv $prefix\.png $prefix\.Non-responder.png
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

