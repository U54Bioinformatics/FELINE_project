#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=9:00:00
#SBATCH --output=betsy_WGS_CNA_Sequenza_FELINE_step2.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


export PATH=$PATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/
export PATH=$PATH:/home/jichen/software/BETSY/install/bin:~/software/BETSY/install/lib64/python2.7/site-packages/genomicode/bin/
export PYTHONPATH=$PYTHONPATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages:/home/jichen/software/Python_lib/lib/python2.7/site-packages/:/home/jichen/software/BETSY/install/lib/python2.7/site-packages
export R_LIBS=/home/jichen/software/BETSY/install/envs/sequenza_2_1_2/lib/R/library
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jichen/software/BETSY/install/envs/sequenza/lib

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
#GTF=/home/jichen/Projects/Database/Reference/hg19/Homo_sapiens.GRCh37.87.chr.gtf


prefix=FELINE
#Segmentation with best fit model
if [ true == true ]; then
#total 1.0 mins on apollo
betsy_run.py --environment coh_slurm --receipt betsy_WGS_CNA_Sequenza_step2a_receipt.txt --num_cores 16 \
    --network_json betsy_WGS_CNA_Sequenza_step2a_network_json.txt --network_png betsy_WGS_CNA_Sequenza_step2a_network.pdf \
    --input SequenzaResults --input_file $prefix\_cnv_sequenza \
    --input ReferenceGenome --input_file $GENOME/Broad.GRCh37/GRCh37/Homo_sapiens_assembly19.fasta \
    --input SequenzaModelSelectionFile --input_file $prefix\_cnv_sequenza/models.highest_probability.txt \
    --output ChromosomalSegmentSignalFile --output_file $prefix\_cnv_sequenza_ChromSegmentSig.txt \
    --dattr ChromosomalSegmentSignalFile.has_multiple_models=no \
    --dattr ChromosomalSegmentSignalFile.evenly_spaced=yes \
    --mattr cn_header=CNt \
    --mattr discard_chrom_with_prefix=GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605,M \
    --run

fi

#get gene copy numbers
if [ true == true ]; then
#total 2.0 mins on apollo
betsy_run.py --environment coh_slurm --receipt betsy_WGS_CNA_Sequenza_step2b_receipt.txt --num_cores 16 \
    --network_json betsy_WGS_CNA_Sequenza_step2b_network_json.txt --network_png betsy_WGS_CNA_Sequenza_step2b_network.pdf \
    --input SequenzaResults --input_file $prefix\_cnv_sequenza \
    --input SequenzaModelSelectionFile --input_file $prefix\_cnv_sequenza/models.highest_probability.txt \
    --input GTFGeneModel --input_file $GTF \
    --output SignalFile --output_file $prefix\_cnv_sequenza_copy_number_gene.txt \
    --dattr SignalFile.has_entrez_gene_info=yes \
    --dattr SignalFile.preprocess=copy_number \
    --mattr facets_gbuild=hg19 \
    --mattr cn_header=CNt \
    --mattr discard_chrom_with_prefix=GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605,M \
    --run
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

