#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=180G
#SBATCH --time=90:00:00
#SBATCH --output=betsy_WGS_CNA_FACETS_FELINE_step1.sh.%A_%a.stdout
#SBATCH -p abild
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

GENOME=/home/jichen/Projects/Database/Genomes
CANCER_GENES=cancer_genes.txt
COSMIC=cosmic.v79.grch37.mutation_data.txt.gz


#Module: compile_pileup_for_facets
#------------------------------------------------------------------------
#BamFolder
#    processing_stage=done
#    sorted=coordinate
#    indexed=yes
#NormalCancerFile
#->
#FACETSPileupFolder

#Attributes:
#facets_snps_vcf     VCF file containing SNPs across the genome, e.g.
#    dbSNP.  Must be sorted.
#facets_min_base_quality(default 20)
#facets_min_map_quality(default 30)
#facets_pseudo_snps  At this number of positions, if there is no SNP,
#    insert a pseudo-SNP. (default 100)
#facets_min_read_counts_germline(default 20)
#facets_min_read_counts_tumor(default 0)
prefix=FELINE
if true; then
betsy_run.py --environment coh_slurm --receipt betsy_WGS_CNA_FACETS_step1a_receipt.txt --num_cores 16 \
    --network_json betsy_WGS_CNA_FACETS_step1a_network_json.txt --network_png betsy_WGS_CNA_FACETS_step1a_network.pdf \
    --input BamFolder --input_file $prefix\_bam \
    --dattr BamFolder.sorted=coordinate \
    --dattr BamFolder.indexed=yes \
    --dattr BamFolder.processing_stage=done \
    --input NormalCancerFile --input_file $prefix\_normal_cancer.txt \
    --output FACETSPileupFolder --output_file $prefix\_FACETSPileupFolder \
    --mattr facets_snps_vcf=$GENOME/Broad.GRCh37/GRCh37/dbsnp_138.b37.vcf.gz \
    --run
fi

#========================================================================
#Module: estimate_copy_number_with_facets
#------------------------------------------------------------------------
#FACETSPileupFolder
#->
#FACETSResults
#
#Attributes:
#facets_gbuild       Genome used for alignment.  Valid values are:
#    hg38, hg19, hg18, mm10, mm9
#facets_critical_valueLower facets_critical_value means higher
#    sensitivity for small changes. (default 150)
#facets_ndepth       Minimum normal sample depth. (default 35)
#facets_nbhd         Window size. (default 250)
#facets_nhet         Minimum number of heterozygous snps in a segment
#    use for bivariate t-statistic. (default 15)
if false; then
betsy_run.py --environment coh_slurm --receipt betsy_WGS_CNA_FACETS_step1b_receipt.txt --num_cores 16 \
    --network_json betsy_WGS_CNA_FACETS_step1b_network_json.txt --network_png betsy_WGS_CNA_FACETS_step1b_network.pdf \
    --input FACETSPileupFolder --input_file $prefix\_FACETSPileupFolder \
    --output FACETSResults --output_file $prefix\_FACETSResults \
    --mattr facets_gbuild=hg19 \
    --run
fi

#if false; then
# parameters setting see this issue
#snp.nbhd is chosen based on insert size. The default is 250 which works well for WES samples with median insert size of 250-300. If your samples have smaller median insert sizes you should set them lower to say 150. Too low a value (relative to insert size) can induce auto correlation in the data which affect the fit.
#ndepth is chosen based is normal coverage. With 100x default 35 will be fine.
#min.nhet is the minimum number of hets in a segment that is used for fitting the allele specific integer copy numbers.
#The key parameter is cval which is used for segmentation. Higher value leads to fewer segments.
#In the end the fit depends on the how well the genome is covered and more importantly the expected number of het snps. The fewer the loci/hets the less robust the coverage. Use the plotSample command to visualize the data and see how noisy it is.
# https://github.com/mskcc/facets/issues/81
#total 38.1 mins on apollo
#fi

round=0
for cval in 100 150 200 300 400 500 600 700 800 900 1000
do
   for ndepth in 35 50
   do
       for nbhd in 150 250
       do
           for nhet in 15
           do

let "round=round+1"
betsy_run.py --environment coh_slurm --receipt betsy_WGS_CNA_FACETS_step1b_receipt.txt --num_cores 16 \
    --network_json betsy_WGS_CNA_FACETS_step1b_network_json.txt --network_png betsy_WGS_CNA_FACETS_step1b_network.pdf \
    --input FACETSPileupFolder --input_file $prefix\_FACETSPileupFolder \
    --output FACETSResults --output_file $prefix\_FACETSResults_parameter$round \
    --mattr facets_gbuild=hg19 \
    --mattr facets_critical_value=$cval \
    --mattr facets_ndepth=$ndepth \
    --mattr facets_nbhd=$nbhd \
    --mattr facets_nhet=$nhet \
    --run
           done
       done
   done           
done

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

