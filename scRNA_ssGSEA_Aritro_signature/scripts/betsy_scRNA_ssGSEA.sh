#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=betsy_scRNA_ssGSEA.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh

#memory and time
#FEL013 use 16 cpu, 1h40min, 42G mem 
#FEL017-27 use 4 cpu, 1h40min-9h, 8-160G mem depending on cell numbers (1k-40k). 

#input count file: use processed file by scRNA_Seurat_00origial_patient_obj.R. This file has only
#one patient. Sometime file "FEL017P109_10x_count01.noaggr.counts.txt" has multiple patients causing
#high memory usage and unwanted results. 


export PATH=$PATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/
#export PATH=$PATH:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages:/home/jichen/software/Python_lib/lib/python2.7/site-packages/:/home/jichen/software/BETSY/install/envs/ssGSEA/lib/python2.7/site-packages
#export PYTHONPATH=$PYTHONPATH:/home/jichen/software/BETSY/install/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/lib/python2.7:/home/jichen/software/BETSY/install/envs/ssGSEA/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/envs/ssGSEA/lib/python2.7
export R_LIBS=/home/jichen/software/BETSY/install/envs/ssGSEA/lib/R/library/

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

#cut into one sequence one file
#perl ~/software/bin/fastaDeal.pl --cuts 1 ~/Projects/Database/Reference/hg19/ucsc.hg19.fasta > hg19.chr1.fasta
#cp ucsc.hg19.fasta.cut/ucsc.hg19.fasta.01 hg19.chrM.fasta
#GENOME=/home/jichen/Projects/Breast/WGS/MDS_Kahn/Reference/hg19
GENOME=/home/jichen/Projects/Database/Reference
Database=/home/jichen/Projects/Database

platform=10x
prefix=`cat patients.list | head -n $N | tail -n 1`
#prefix=`cat patients_chunks.list | head -n $N | tail -n 1`
#prefix=$sample\_$platform
#total 10.1 hrs on apollo
betsy_run.py --environment coh_slurm --receipt betsy_scRNA_ssGSEA_receipt.txt --num_cores $CPU \
    --network_json betsy_scRNA_ssGSEA_network_json.txt --network_png betsy_scRNA_ssGSEA_network.pdf \
    --input SignalFile --input_file $prefix\.zinbwave.normalized.txt \
    --dattr SignalFile.preprocess=tpm \
    --dattr SignalFile.has_human_gene_info=no\
    --dattr SignalFile.logged=not_applicable \
    --output ssGSEAResults --output_file $prefix\.ssGSEA \
    --mattr geneset_file=AN_GeneSets.gmt \
    --mattr gsva_min_geneset_size=10 \
    --run
 
end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

