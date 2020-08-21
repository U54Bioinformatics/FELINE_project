#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=80G
#SBATCH --time=60:00:00
#SBATCH --output=betsy_WGS_variants_FELINE.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


export PATH=$PATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/
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

#grep "Wed" betsy_WGS_variants_FELINE.sh.1782001_2.stdout | cut -d" " -f5 | awk '{print "rm -R "$1"_*"}' > clean_WES_cache.sh

#cut into one sequence one file
#perl ~/software/bin/fastaDeal.pl --cuts 1 ~/Projects/Database/Reference/hg19/ucsc.hg19.fasta > hg19.chr1.fasta
#cp ucsc.hg19.fasta.cut/ucsc.hg19.fasta.01 hg19.chrM.fasta
#GENOME=/home/jichen/Projects/Breast/WGS/MDS_Kahn/Reference/hg19
GENOME=/home/jichen/Projects/Database/Genomes
CANCER_GENES=cancer_genes.txt
COSMIC=cosmic.v79.grch37.mutation_data.txt.gz

#FILE=`cat FELINE_patients.14_23_26_27_37_39.list | head -n $N | tail -n 1`
#FILE=`cat FELINE_patients.44_45_46.list | head -n $N | tail -n 1`
#FILE=`cat FELINE_patients.25.list | head -n $N | tail -n 1`
#FILE=`cat FELINE_patients.25_31_43.list | head -n $N | tail -n 1`
FILE=`cat FELINE_patients.list | head -n $N | tail -n 1`

prefix=$FILE
echo $prefix
#total 1 day, 5.7 hours on apollo
#Originalbam files should be in the same folder because BESTY trace back to the folder and use
#folder name from one sample for all others
betsy_run.py --environment coh_slurm --receipt betsy_WGS_variants_receipt.txt --num_cores $CPU \
    --network_json betsy_WGS_variants_network_json.txt --network_png betsy_WGS_variants_network.pdf \
    --input BamFolder --input_file $prefix\_bam \
    --dattr BamFolder.sorted=coordinate \
    --dattr BamFolder.has_read_groups=yes \
    --dattr BamFolder.duplicates_marked=yes \
    --dattr BamFolder.indel_realigned=yes \
    --dattr BamFolder.base_quality_recalibrated=yes \
    --dattr BamFolder.indexed=yes \
    --dattr BamFolder.aligner=bwa_mem \
    --input SampleGroupFile --input_file $prefix\_sample.txt \
    --input ReferenceGenome --input_file $GENOME/Broad.GRCh37/GRCh37/Homo_sapiens_assembly19.fasta \
    --input NormalCancerFile --input_file $prefix\_normal_cancer.txt \
    --output SimpleVariantMatrix --output_file $prefix\_variant_calls \
    --also_save_highest ManyCallerVCFFolders,$prefix\_call05 \
    --dattr SimpleVariantMatrix.duplicates_marked=yes \
    --dattr SimpleVariantMatrix.caller_suite=lance3 \
    --mattr wgs_or_wes=wes \
    --mattr mutect_dbsnp_vcf=$GENOME/Broad.GRCh37/GRCh37/dbsnp_138.b37.vcf.gz \
    --mattr mutect_cosmic_vcf=$GENOME/COSMIC/b37_cosmic_v54_120711.vcf \
    --dattr SimpleVariantMatrix.filtered_calls=yes \
    --mattr filter_by_min_total_reads=10 \
    --mattr filter_by_min_alt_reads=2 \
    --mattr filter_by_min_vaf=0.01 \
    --mattr filter_by_min_callers=1 \
    --dattr SimpleVariantMatrix.annotated_with_annovar=yes \
    --mattr annovar_buildver=hg19 \
    --dattr SimpleVariantMatrix.annotated_with_snpeff=yes \
    --mattr snpeff_genome=hg19 \
    --dattr SimpleVariantMatrix.variant_annotations=cancer_cosmic \
    --mattr cosmic_variants_file=$GENOME/008_Cancer_Genes/${COSMIC} \
    --mattr cancer_genes_file=$GENOME/008_Cancer_Genes/${CANCER_GENES} \
    --dattr SimpleVariantMatrix.with_coverage=yes \
    --mattr discard_chrom_with_prefix=GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605,M \
    --run 

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

