#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=80G
#SBATCH --time=40:00:00
#SBATCH --output=betsy_WGS_CNA_Sequenza_FELINE_step3.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


export PATH=$PATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/
export PATH=$PATH:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages:/home/jichen/software/Python_lib/lib/python2.7/site-packages/:/home/jichen/software/BETSY/install/lib/python2.7/site-packages
export R_LIBS=/home/jichen/software/BETSY/install/envs/CNA/lib/R/library

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
GENOME=/home/jichen/Projects/Database/Genomes
CANCER_GENES=cancer_genes.txt
COSMIC=cosmic.v79.grch37.mutation_data.txt.gz
GTF=/home/jichen/Projects/Database/Genomes/annotations/Homo_sapiens.GRCh37.87.gtf

prefix=FELINE

#Gene level copy number
if true; then
#total < 2mins on apollo
betsy_run.py --environment coh_slurm --receipt betsy_WGS_CNA_Sequenza_step3b_receipt.txt --num_cores 16 \
    --network_json betsy_WGS_CNA_Sequenza_step3b_network_json.txt --network_png betsy_WGS_CNA_Sequenza_step3b_network.pdf \
    --input SequenzaResults --input_file $prefix\_cnv_sequenza \
    --input SequenzaModelSelectionFile --input_file $prefix\_cnv_sequenza_model_selection.txt \
    --input ReferenceGenome --input_file $GENOME/GRCh37/genome.fa \
    --input GTFGeneModel --input_file $GTF \
    --output CopyNumberAnalysis --output_file $prefix\_cnv_sequenza_copy_number_analysis.txt \
    --dattr FACETSModelSelectionFile.model_selection=adhoc \
    --mattr facets_gbuild=hg19 \
    --mattr cn_header=CNt \
    --mattr total_cn_header=CNt \
    --mattr minor_cn_header=B \
    --mattr discard_chrom_with_prefix=MT \
    --run
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

