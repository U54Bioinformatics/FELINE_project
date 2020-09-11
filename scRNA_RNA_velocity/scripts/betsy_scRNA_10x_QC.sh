#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=180G
#SBATCH --time=60:00:00
#SBATCH --output=betsy_scRNA_10x_QC.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


#export PATH=$PATH:/home/jichen/software/BETSY/install/envs/scRNA/bin/:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
#export PYTHONPATH=$PYTHONPATH:/home/jichen/software/BETSY/install/lib/python2.7/site-packages:/home/jichen/software/BETSY/install/lib/python2.7:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
#export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library

export PATH=$PATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages:/home/jichen/software/Python_lib/lib/python2.7/site-packages/:/home/jichen/software/BETSY/install/envs/scRNA/lib/python3.7/site-packages
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library:/home/jichen/software/BETSY/install/envs/scRNAcnv/lib/R/library/

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

prefix=`cat patients.list | head -n $N | tail -n 1`
#preprocess 10x data, total 14.8 hrs on apollo
if true ; then
echo "preprocess 10x data ......"
betsy_run.py --environment coh_slurm --receipt betsy_scRNA_10x_QC_receipt01.txt --num_cores $CPU \
    --network_json betsy_scRNA_10x_QC_network_json.txt --network_png betsy_scRNA_10x_QC_network.pdf \
    --input TenXFastqFolder --input_file $prefix\_fastq \
    --input SampleGroupFile --input_file $prefix\_sample.txt \
    --output SignalFile --output_file $prefix\_count01.noaggr.txt \
    --dattr SignalFile.gene_expression_estimator=featurecounts \
    --dattr SignalFile.aligner=star \
    --dattr SignalFile.has_entrez_gene_info=yes \
    --dattr TenXSignalFile.aggregated=no \
    --also_save_highest TenXCountResults,$prefix\_count03 \
    --mattr tenx_refdata_path=$Database/Genomes/refdata-cellranger-GRCh38-3.0.0 \
    --run
fi

#Do the QC on the analysis.
#betsy_run.py --environment coh_slurm --receipt betsy_scRNA_10x_QC_receipt02.txt --num_cores 16 \
#    --input TenXCountResults --input_file $prefix\_count03.txt \
#    --output TenXQCFolder --output_file $prefix\_qc01 \
#    --dattr TenXQCFolder.aggregated=no


#QC on 10x data, total 6.4 hrs on apollo
if [ true == true ]; then
echo "QC on 10x data ......"
betsy_run.py --environment coh_slurm --receipt betsy_scRNA_10x_QC_receipt02.txt --num_cores $CPU \
     --input SignalFile --input_file $prefix\_count01.noaggr.txt \
     --dattr SignalFile.preprocess=counts \
     --dattr SignalFile.logged=no \
     --dattr SignalFile.gene_expression_estimator=featurecounts \
     --output scRNAMetadata --output_file $prefix\_meta \
     --mattr sample_names_are_10x_format=yes \
     --run

xls2txt $prefix\_meta/summary.xls > $prefix\_meta/summary.txt

fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

