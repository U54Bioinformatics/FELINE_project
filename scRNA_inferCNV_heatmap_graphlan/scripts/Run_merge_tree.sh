#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --output=Run_merge_tree.sh.%A_%a.stdout
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

echo "sort patients"
# select patient with cancer cells pre and post treatment
/home/jichen/software/BETSY/install/envs/R_combine_pdf_figures/bin/Rscript Run_merge_tree_helper.R

echo "prepare script"
prefix=FELINE_scRNA_tree_24patient

# move individual tree figure into MEDALT_tree_png
if [ ! -e MEDALT_tree_png ]; then
    mkdir MEDALT_tree_png
    mv *.tree.png MEDALT_tree_png
fi

# Non-responder
if true; then
response=Non-responder
python Run_merge_tree_helper.py --input FELINE_patient_info.sorted.txt --project $prefix --response $response
/home/jichen/software/BETSY/install/envs/R_combine_pdf_figures/bin/Rscript $prefix\.$response\.R
convert -density 100 -quality 100 $prefix\.$response\.pdf $prefix\.$response\.png
fi

# Responder
if true; then
response=Responder
python Run_merge_tree_helper.py --input FELINE_patient_info.sorted.txt --project $prefix --response $response
/home/jichen/software/BETSY/install/envs/R_combine_pdf_figures/bin/Rscript $prefix\.$response\.R
convert -density 100 -quality 100 $prefix\.$response\.pdf $prefix\.$response\.png
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

