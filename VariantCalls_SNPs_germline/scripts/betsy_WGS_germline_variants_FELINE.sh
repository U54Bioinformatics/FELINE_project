#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=80G
#SBATCH --time=40:00:00
#SBATCH --output=betsy_WGS_germline_variants_FELINE.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


export PATH=/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/:$PATH
export PYTHONPATH=/home/jichen/software/Python_lib64/lib64/python2.7/site-packages:/home/jichen/software/Python_lib/lib/python2.7/site-packages/:/home/jichen/software/BETSY/install/lib/python2.7/site-packages
module load singularity/3.5.3


#export R_LIBS=/home/hmirsafian/miniconda3/envs/mypython/lib/R/library
#export PATH=/home/hmirsafian/Software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin:$PATH
#export PYTHONPATH=/home/hmirsafian/Software/Python_lib64/lib64/python2.7/site-packages/genomicode:/home/hmirsafian/Software/Python_lib/lib/python2.7/site-packages:/home/hmirsafian/miniconda3/envs/mypython/lib/python2.7/site-packages
#export LD_LIBRARY_PATH=/home/hmirsafian/miniconda3/envs/mypython/lib/:/home/hmirsafian/miniconda3/envs/mypython/lib/R/lib


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
bed=/home/jichen/Projects/Breast/scRNA/Data/Agilent_SureSelect_XT_Human_v7.nochr.bed

FILE=`cat FELINE_patients.4.list | head -n $N | tail -n 1`

prefix=$FILE
echo $prefix
#total 1 day, 5.7 hours on apollo
betsy_run.py --environment coh_slurm --receipt betsy_WGS_germline_variants_receipt.txt --num_cores $CPU \
    --network_json betsy_WGS_germline_variants_network_json.txt --network_png betsy_WGS_germline_variants_network.pdf \
    --input BamFolder --input_file $prefix\_bam \
    --dattr BamFolder.each_file_contains=one_sample \
    --dattr BamFolder.sorted=coordinate \
    --dattr BamFolder.has_read_groups=yes \
    --dattr BamFolder.duplicates_marked=yes \
    --dattr BamFolder.indel_realigned=yes \
    --dattr BamFolder.base_quality_recalibrated=yes \
    --dattr BamFolder.indexed=yes \
    --dattr BamFolder.aligner=bwa_mem \
    --dattr BamFolder.adapters_trimmed=yes \
    --input SampleGroupFile --input_file $prefix\_sample.txt \
    --input ReferenceGenome --input_file $GENOME/Broad.GRCh37/GRCh37/Homo_sapiens_assembly19.fasta \
    --output VariantCallingAnalysis --output_file $prefix\_germline_calls \
    --dattr VariantCallingAnalysis.duplicates_marked=yes \
    --dattr VariantCallingAnalysis.caller_suite=general \
    --dattr VCFFolder.aligner=bwa_mem \
    --mattr wgs_or_wes=wes \
    --dattr VariantCallingAnalysis.filtered_calls=yes \
    --mattr filter_by_min_total_reads=10 \
    --mattr filter_by_min_alt_reads=2 \
    --mattr filter_by_min_vaf=0.05 \
    --mattr filter_by_min_callers=1 \
    --dattr VariantCallingAnalysis.annotated_with_annovar=yes \
    --mattr annovar_buildver=hg19 \
    --dattr VariantCallingAnalysis.annotated_with_snpeff=yes \
    --mattr snpeff_genome=GRCh37.75 \
    --dattr VariantCallingAnalysis.variant_annotations=cancer_cosmic \
    --mattr cancer_genes_file=$GENOME/008_Cancer_Genes/${CANCER_GENES} \
    --mattr cosmic_variants_file=$GENOME/008_Cancer_Genes/${COSMIC} \
    --dattr VariantCallingAnalysis.with_coverage=yes \
    --mattr tmb_targeted_regions=$bed \
    --mattr discard_chrom_with_prefix=GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605,M \
    --run

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

