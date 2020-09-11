#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --output=run_vaf2maf.sh.%A_%a.stdout
#SBATCH -p fast,all
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/lib/R/library:$R_LIBS
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/lib/python2.7/site-packages/$PYTHONPATH
export PERL5LIB=$PERL5LIB:~/software/BETSY/install/envs/DNA_analysis_WES/lib/site_perl/5.26.2/
export PATH=/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/bin:$PATH
export VEP_DATA=/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/VEP_lib/hg19/vep/homo_sapiens/86_GRCh37


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

vep_bin=/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/bin/
vep_cache=/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/VEP_lib/hg19/vep/
#filter_vcf=$vep_cache/homo_sapiens/86_GRCh37/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
filter_vcf=0

#vcf2maf.pl --input-vcf tests/test.vcf --output-maf tests/test.vcf.maf --vep-path $vep_bin --vep-data $vep_cache --ref-fasta $vep_cache/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --filter-vcf $vep_cache/homo_sapiens/86_GRCh37/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz 
maf2maf.pl --input-maf FELINE.vep.tsv --output-maf FELINE.vep.maf --vep-path $vep_bin --vep-data $vep_cache --ref-fasta $vep_cache/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --filter-vcf $filter_vcf

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

