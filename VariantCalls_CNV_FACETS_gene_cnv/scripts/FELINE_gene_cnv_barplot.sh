#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=4:00:00
#SBATCH --output=FELINE_gene_cnv_barplot.sh.%A_%a.stdout
#SBATCH -p abild,all
#SBATCH --workdir=./

export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library/

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
echo "Command line N: $1"

# summarize by per sample/biopsy
#/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript FELINE_gene_cnv_barplot.R
# summarize by per patient
/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript FELINE_gene_cnv_barplot.patient.R
# do fisher exact test for gene cnv between responder vs. nonresponder
# output 1 is FELINE_gene_cnv_barplot.patient.STable_fisher_test.txt
# output 2 is FELINE_gene_cnv_barplot.patient.STable_allCNV.txt, a patient gene cnv calls 
/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript FELINE_gene_cnv_barplot.patient.fisher_test.R 
# add patient three version of IDs to cnv calls for Jason
# output is FELINE_gene_cnv_barplot.patient.STable_allCNV.three_IDs.txt
/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript add_patient_id.R




end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

