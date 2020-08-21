#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=20:00:00
#SBATCH --output=betsy_WGS_pyclone.sh.%A_%a.stdout
#SBATCH -p abild,all
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


export PATH=/home/jichen/software/Python_20200309lib64/lib64/python2.7/site-packages/genomicode/bin/:$PATH
export PYTHONPATH=~/software/BETSY/install/envs/pyclone/lib/python2.7/site-packages/:/home/jichen/software/Python_20200309lib64/lib64/python2.7/site-packages:/home/jichen/software/Python_20200309lib/lib/python2.7/site-packages/:/home/jichen/software/BETSY/install/lib/python2.7/site-packages
export R_LIBS=/home/jichen/software/BETSY/install/lib/R/library


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
#GTF=/home/jichen/Projects/Database/Reference/hg19/Homo_sapiens.GRCh37.87.chr.gtf
GTF=/home/jichen/Projects/Database/Reference/GTF_analysis/Homo_sapiens.GRCh37.87.chr.gtf
#GTF_nochr=/home/jichen/Projects/Database/Reference/GTF_analysis/Homo_sapiens.GRCh37.87.gtf
GTF_nochr=/home/jichen/Projects/Database/Genomes/annotations/Homo_sapiens.GRCh37.87.gtf

prefix=`cat FELINE_patients.list | head -n $N | tail -n 1`
#prefix=`cat FELINE_patients.drop.list | head -n $N | tail -n 1`
#prefix=`cat FELINE_patients.14_23_26_27_37_39.list | head -n $N | tail -n 1`
#prefix=FELINE

# run option
pyclone_sequenza_pre=false
pyclone_sequenza_run=false
pyclone_FACETS_pre=false
pyclone_FACETS_run=true

#filter variants SNP/INDEL
min_total_read=20
min_alt_reads=5
min_vaf=0.05
min_callers=1 
max_copy=none # none means do not use this seting. use all cnv regions
cn=no
germline=$prefix\_G
pyclone_iterations=10000 # 1000 for testing run and 10000 for final run

echo "Pamameters: min_total_read, min_alt_reads, min_vaf, min_callers, pyclone_iterations"
echo "$min_total_read, $min_alt_reads, $min_vaf, $min_callers, $pyclone_iterations"
echo "Pamameters: pyclone_sequenza_pre, pyclone_sequenza_run, pyclone_FACETS_pre, pyclone_FACETS_run"
echo "$pyclone_sequenza_pre, $pyclone_sequenza_run, $pyclone_FACETS_pre, $pyclone_FACETS_run"

if [ ! -e $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.txt ] ; then
# remove any germline variants in somatic mutations
#slice_svm.py --keep_samples_startswith $germline $prefix\_variant_calls > $prefix\_variant_calls_Germline
#slice_svm.py --filter_min_alt_reads 1 $prefix\_variant_calls_Germline > $prefix\_variant_calls_Germline_variant
#slice_svm.py --convert_to_position_file $prefix\_variant_calls_Germline_variant | grep -v "chrom" | awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4}' > $prefix\_variant_calls_Germline_variant.bed
    if [ $prefix == 'FEL044' ] || [ $prefix == 'FEL045' ] || [ $prefix == 'FEL046' ] || [ $prefix == 'FEL036' ] ; then
        python Prepare_somatic_variant_Germline_variants_FELINE.py --input $prefix\_variant_calls --germline_col 15
    else
        python Prepare_somatic_variant_Germline_variants_FELINE.py --input $prefix\_variant_calls --germline_col 16
    fi
script=~/software/BETSY/scripts/
$script/slice_svm.py --discard_by_coord $prefix\_variant_calls_Germline_variant.bed $prefix\_variant_calls > $prefix\_variant_calls_Germline_free
# common dbsnp info
python Prepare_somatic_variant_common_variant_FELINE.py --input $prefix\_variant_calls_Germline_free --dbsnp_col 26
# filter somatic mutations by read depth and var and number of callers
$script/slice_svm.py --remove_sample $germline $prefix\_variant_calls_Germline_free > $prefix\_variant_calls_01_purity40.txt
$script/slice_svm.py --filter_min_total_reads $min_total_read $prefix\_variant_calls_01_purity40.txt > $prefix\_variant_calls_02_purity40_totalRD$min_total_read\.txt
$script/slice_svm.py --filter_min_alt_reads $min_alt_reads $prefix\_variant_calls_02_purity40_totalRD$min_total_read\.txt > $prefix\_variant_calls_03_purity40_totalRD$min_total_read\_minread$min_alt_reads\.txt
$script/slice_svm.py --discard_min_vaf $min_vaf $prefix\_variant_calls_03_purity40_totalRD$min_total_read\_minread$min_alt_reads\.txt > $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.txt
$script/slice_svm.py --filter_min_callers $min_callers $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.txt > $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.txt
$script/slice_svm.py --convert_to_position_file $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.txt | grep -v "chrom" | awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4}' > $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.bed
$script/slice_svm.py --convert_to_position_file $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.txt | grep -v "chrom" | awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4}' > $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.bed
# common dbsnp info
python Prepare_somatic_variant_common_variant_FELINE.py --input $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.txt --dbsnp_col 25
# exonic only
$script/slice_svm.py --exonic_only $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.txt > $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.exonic.txt
$script/slice_svm.py --exonic_only $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.txt > $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.exonic.txt
python Prepare_somatic_variant_VCF_FELINE.py --input $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.exonic.txt --vcf_header vcf_header.vcf
python Prepare_somatic_variant_VCF_FELINE.py --input $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.txt --vcf_header vcf_header.vcf
python Prepare_somatic_variant_VCF_FELINE.py --input $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.txt --vcf_header vcf_header.vcf
# common dbsnp info
python Prepare_somatic_variant_common_variant_FELINE.py --input $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.exonic.txt --dbsnp_col 25
# For FACETS
# minVAF
sed 's/chr//g' $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.txt > $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.FACETS.txt
sed 's/chr//g' $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.exonic.txt > $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.exonic.FACETS.txt
python Prepare_somatic_variant_VCF_FELINE.py --input $prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.FACETS.txt --vcf_header vcf_header.vcf
# mincaller
sed 's/chr//g' $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.txt > $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.FACETS.txt
sed 's/chr//g' $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.exonic.txt > $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.exonic.FACETS.txt
python Prepare_somatic_variant_VCF_FELINE.py --input $prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.FACETS.txt --vcf_header vcf_header.vcf
fi

# when 1 will not use mincaller filter for pyclone
if [ $min_callers == 1 ]; then
    echo "filter min_total_read, min_alt_reads, and min_vaf"
    # filter min_total_read, min_alt_reads, and min_vaf
    input_var=$prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.txt
    input_var_FACETS=$prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.FACETS.txt
    if [ $prefix == 'FEL045' ] || [ $prefix == 'FEL045' ]; then
        input_var=$prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.exonic.txt
        input_var_FACETS=$prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.exonic.FACETS.txt
    fi
fi

# when 2 will use mincaller filter for pyclone
if [ $min_callers == 2 ]; then
    echo "filter more with min_callers"
    # filter more with min_callers
    input_var=$prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.txt
    input_var_FACETS=$prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.FACETS.txt
    if [ $prefix == 'FEL045' ] || [ $prefix == 'FEL045' ]; then
        input_var=$prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.exonic.txt
        input_var_FACETS=$prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.exonic.FACETS.txt
    fi
fi

if [ ! -e $prefix\_pyclone_mutation_sequenza ] && $pyclone_sequenza_pre ; then
#analyze mutation numbers/quality for pyclone, sequenza (total 5 secs on apollo)
#input_var=$prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.txt
#input_var=$prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.txt
betsy_run.py --environment coh_slurm --receipt betsy_WGS_pyclone_step1a_receipt.txt --num_cores $CPU \
    --network_json betsy_WGS_pyclone_step1a_network_json.txt --network_png betsy_WGS_pyclone_step1a_network.pdf \
    --input SimpleVariantMatrix --input_file $input_var \
    --dattr SimpleVariantMatrix.with_coverage=yes \
    --input SequenzaResults --input_file $prefix\_cnv_sequenza \
    --input SequenzaModelSelectionFile --input_file $prefix\_cnv_sequenza/models.highest_probability.txt \
    --output PyCloneMutationsFolder --output_file $prefix\_pyclone_mutation_sequenza \
    --mattr cn_header=CNt \
    --mattr total_cn_header=CNt \
    --mattr minor_cn_header=B \
    --mattr min_coverage_for_pyclone_variants=$min_total_read \
    --mattr use_only_consistent_cn=$cn\
    --mattr pyclone_iterations=$pyclone_iterations \
    --run
    xls2txt $prefix\_pyclone_mutation_sequenza/num_mutations.xls > $prefix\_pyclone_mutation_sequenza/num_mutations.txt
fi
#do not see max_copynum will use all cnv regions
#--mattr max_copynum_for_pyclone=$max_copy \

if [ ! -e $prefix\_pyclone_mutation_FACETS  ] && $pyclone_FACETS_pre ; then
#analyze mutation numbers/quality for pyclone, FACETS (total 5 secs on apollo)
#input_var_FACETS=$prefix\_variant_calls_04_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\.FACETS.txt
#input_var_FACETS=$prefix\_variant_calls_05_purity40_totalRD$min_total_read\_minread$min_alt_reads\_minVAF$min_vaf\_mincaller$min_callers\.FACETS.txt
betsy_run.py --environment coh_slurm --receipt betsy_WGS_pyclone_step1b_receipt.txt --num_cores $CPU \
    --network_json betsy_WGS_pyclone_step1b_network_json.txt --network_png betsy_WGS_pyclone_step1b_network.pdf \
    --input SimpleVariantMatrix --input_file $input_var_FACETS \
    --dattr SimpleVariantMatrix.with_coverage=yes \
    --input FACETSResults --input_file $prefix\_FACETSResults \
    --input FACETSModelSelectionFile --input_file $prefix\_FACETSResults/model_selection.txt \
    --output PyCloneMutationsFolder --output_file $prefix\_pyclone_mutation_FACETS \
    --mattr facets_gbuild=hg19 \
    --mattr cn_header=tcn.em \
    --mattr total_cn_header=tcn.em \
    --mattr minor_cn_header=lcn.em \
    --mattr min_coverage_for_pyclone_variants=$min_total_read \
    --mattr use_only_consistent_cn=$cn\
    --mattr pyclone_iterations=$pyclone_iterations \
    --run
    xls2txt $prefix\_pyclone_mutation_FACETS/num_mutations.xls > $prefix\_pyclone_mutation_FACETS/num_mutations.txt
fi



if $pyclone_sequenza_run ; then
#run pyclone analysis, sequenza
betsy_run.py --environment coh_slurm --receipt betsy_WGS_pyclone_step1c_receipt.txt --num_cores $CPU \
    --network_json betsy_WGS_pyclone_step1c_network_json.txt --network_png betsy_WGS_pyclone_step1c_network.pdf \
    --input SimpleVariantMatrix --input_file $input_var \
    --dattr SimpleVariantMatrix.with_coverage=yes \
    --input SequenzaResults --input_file $prefix\_cnv_sequenza \
    --input SequenzaModelSelectionFile --input_file $prefix\_cnv_sequenza/models.highest_probability.txt \
    --input GTFGeneModel --input_file $GTF_nochr \
    --output PyCloneAnalysis --output_file $prefix\_pyclone_analysis_sequenza_$min_total_read\_$min_alt_reads\_$min_vaf\_$min_callers\_$pyclone_iterations\_$max_copy \
    --mattr cn_header=CNt \
    --mattr total_cn_header=CNt \
    --mattr minor_cn_header=B \
    --mattr min_coverage_for_pyclone_variants=$min_total_read \
    --mattr use_only_consistent_cn=$cn\
    --mattr pyclone_iterations=$pyclone_iterations \
    --run

    echo "convert png to pdf"
    pyclone_outdir=$prefix\_pyclone_analysis_sequenza_$min_total_read\_$min_alt_reads\_$min_vaf\_$min_callers\_$pyclone_iterations\_$max_copy
    for fig in $pyclone_outdir/cancer_cell_frequency.hm $pyclone_outdir/cluster.frequency $pyclone_outdir/cluster.num_mutations $pyclone_outdir/pairwise.frequency.hm $pyclone_outdir/robustness.score $pyclone_outdir/variant_allele_frequency.hm
    do
       /home/jichen/software/BETSY/install/bin/convert -density 300 -quality 100 $fig\.png $fig\.pdf
    done

    echo "create mutation table with cellular prevalance (CP)"
    for tp in "cluster" "loci" "old_style"
    do
        echo $tp
        burnin_n=`echo $(($pyclone_iterations/10))`
        thin_n=1
        yaml=$pyclone_outdir/pyclone.output/config.yaml
        table=$pyclone_outdir/tables
        mkdir $table
        echo "$burnin_n, $thin_n, $yaml"
        PyClone build_table --config_file $yaml --out_file $table/pyclone.table_$tp\.txt --table_type $tp --burnin $burnin_n --thin $thin_n
    done

    echo "create final cluster file for clonevol"
    if [ ! -e $pyclone_outdir/tables/pyclone.table_cluster.final.txt ]; then
       python Prepare_somatic_variant_pyclone_cluster.py --input $pyclone_outdir/tables/pyclone.table_cluster.txt
    fi
fi

if $pyclone_FACETS_run ; then
#run pyclone analysis, FACETS
betsy_run.py --environment coh_slurm --receipt betsy_WGS_pyclone_step1d_receipt.txt --num_cores $CPU \
    --network_json betsy_WGS_pyclone_step1d_network_json.txt --network_png betsy_WGS_pyclone_step1d_network.pdf \
    --input SimpleVariantMatrix --input_file $input_var_FACETS \
    --dattr SimpleVariantMatrix.with_coverage=yes \
    --input FACETSResults --input_file $prefix\_FACETSResults \
    --input FACETSModelSelectionFile --input_file $prefix\_FACETSResults/model_selection.txt \
    --input GTFGeneModel --input_file $GTF_nochr \
    --output PyCloneAnalysis --output_file $prefix\_pyclone_analysis_FACETS_$min_total_read\_$min_alt_reads\_$min_vaf\_$min_callers\_$pyclone_iterations\_$max_copy \
    --mattr facets_gbuild=hg19 \
    --mattr cn_header=tcn.em \
    --mattr total_cn_header=tcn.em \
    --mattr minor_cn_header=lcn.em \
    --mattr min_coverage_for_pyclone_variants=$min_total_read \
    --mattr use_only_consistent_cn=$cn\
    --mattr pyclone_iterations=$pyclone_iterations \
    --run

    echo "convert png to pdf"
    pyclone_outdir=$prefix\_pyclone_analysis_FACETS_$min_total_read\_$min_alt_reads\_$min_vaf\_$min_callers\_$pyclone_iterations\_$max_copy
    for fig in $pyclone_outdir/cancer_cell_frequency.hm $pyclone_outdir/cluster.frequency $pyclone_outdir/cluster.num_mutations $pyclone_outdir/pairwise.frequency.hm $pyclone_outdir/robustness.score $pyclone_outdir/variant_allele_frequency.hm
    do
       /home/jichen/software/BETSY/install/bin/convert -density 300 -quality 100 $fig\.png $fig\.pdf
    done

    echo "create mutation table with cellular prevalance (CP)"
    for tp in "cluster" "loci" "old_style"
    do
        echo $tp
        burnin_n=`echo $(($pyclone_iterations/10))`
        thin_n=1
        yaml=$pyclone_outdir/pyclone.output/config.yaml
        table=$pyclone_outdir/tables
        mkdir $table
        echo "$burnin_n, $thin_n, $yaml"
        PyClone build_table --config_file $yaml --out_file $table/pyclone.table_$tp\.txt --table_type $tp --burnin $burnin_n --thin $thin_n
    done

    echo "create final cluster file for clonevol"
    if [ ! -e $pyclone_outdir/tables/pyclone.table_cluster.final.txt ]; then
       python Prepare_somatic_variant_pyclone_cluster.py --input $pyclone_outdir/tables/pyclone.table_cluster.txt
    fi
fi


#    --mattr min_coverage_for_pyclone_variants=60 \
#    --mattr max_copynum_for_pyclone=4 \
#    --mattr use_only_consistent_cn=yes \

#interpre pyclone subclones
#: <<'END'
#betsy_run.py --environment coh_slurm --receipt betsy_WGS_pyclone_step1c_receipt.txt --num_cores $CPU \
#    --network_json betsy_WGS_pyclone_step1c_network_json.txt --network_png betsy_WGS_pyclone_step1c_network.pdf \
#    --input PyCloneAnalysis --input_file $prefix\_pyclone_analysis \
#    --output PyCloneSubclones --output_file $prefix\_pyclone_subclones \
#    --run
#    --mattr min_coverage_for_pyclone_variants=60 \
#    --mattr max_copynum_for_pyclone=4 \
#    --mattr use_only_consistent_cn=yes \
#END



end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

