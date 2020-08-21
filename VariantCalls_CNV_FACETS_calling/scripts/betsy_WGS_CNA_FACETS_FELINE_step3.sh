#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=80G
#SBATCH --time=9:00:00
#SBATCH --output=betsy_WGS_CNA_FACETS_FELINE_step3.sh.%A_%a.stdout
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

GENOME=/home/jichen/Projects/Database/Genomes
CANCER_GENES=cancer_genes.txt
COSMIC=cosmic.v79.grch37.mutation_data.txt.gz
GTF=/home/jichen/Projects/Database/Genomes/annotations/Homo_sapiens.GRCh37.87.gtf

prefix=FELINE

###############################################
#1. FELINE_FACETSResults_parameter0/1/2/3/4 are results with different parameters
#2. FELINE_FACETSResults (results with different parameters in this folder): facets.cval_100.ndepth_35.nbhd_250.nhet_15, facets.cval_150.ndepth_35.nbhd_250.nhet_15
#3. FELINE_FACETSResults/summary.txt with all parameters
#Sample  Purity  Ploidy  Diploid LogR    cval    ndepth  snp.nbhd        min.nhet
#FEL044_M        0.720670815287  4.04488576155   -0.796468798153 150     35      250     15
#FEL044_S        0.659488422167  4.00590010471   -0.73242893796  150     35      250     15
#FEL045_M        0.431174985364  2.01476352294   -0.00458456226591       150     35      250     15
#FEL045_S        0.386369957583  2.02285408001   -0.00635557249493       150     35      250     15
#FEL046_M        0.685027907838  2.11442736362   -0.055463574561 150     35      250     15
#FEL046_S                2.0     0.0532762568219 150     35      250     15
#FEL044_M        0.722360518142  4.04028455626   -0.796523420205 100     35      250     15
#FEL044_S        0.689025546264  4.0559643852    -0.772566413671 100     35      250     15
#FEL045_M        0.432532880086  2.01461277035   -0.00455209053921       100     35      250     15
#FEL045_S        0.386149880936  2.02315413377   -0.00643517018722       100     35      250     15
#FEL046_M        0.687263024295  2.11418639921   -0.0555261518779        100     35      250     15
#FEL046_S                2.0     0.0511863373233 100     35      250     15
#4. FELINE_FACETS_model_selection.txt, edit this file for best model and rerun segment and gene copy number
###############################################
run=false # set to false, so we can use a manual curated best model selection
data=runmerge
#Segment all
if $run; then
#total < 2mins on apollo
betsy_run.py --environment coh_slurm --receipt betsy_WGS_CNA_FACETS_step3a_receipt.txt --num_cores 16 \
    --network_json betsy_WGS_CNA_FACETS_step3a_network_json.txt --network_png betsy_WGS_CNA_FACETS_step3a_network.pdf \
    --input FACETSResults --input_file $prefix\_FACETSResults_$data \
    --input ReferenceGenome --input_file $GENOME/GRCh37/genome.fa \
    --output ChromosomalSegmentSignalFile --output_file $prefix\_FACETS_$data\_ChromSegmentSig.cor.txt \
    --also_save_lowest FACETSModelSelectionFile,$prefix\_FACETS_$data\_model_selection.txt \
    --dattr ChromosomalSegmentSignalFile.has_multiple_models=no \
    --dattr ChromosomalSegmentSignalFile.split_chromosomes=yes \
    --dattr ChromosomalSegmentSignalFile.evenly_spaced=yes \
    --dattr FACETSModelSelectionFile.model_selection=best_correlation \
    --mattr facets_gbuild=hg19 \
    --mattr cn_header=tcn.em \
    --mattr discard_chrom_with_prefix=MT \
    --run
fi

#Gene level copy number
if true; then
#total < 2mins on apollo
betsy_run.py --environment coh_slurm --receipt betsy_WGS_CNA_FACETS_step3b_receipt.txt --num_cores 16 \
    --network_json betsy_WGS_CNA_FACETS_step3b_network_json.txt --network_png betsy_WGS_CNA_FACETS_step3b_network.pdf \
    --input FACETSResults --input_file $prefix\_FACETSResults_$data \
    --input FACETSModelSelectionFile --input_file $prefix\_FACETS_$data\_model_selection.txt \
    --input ReferenceGenome --input_file $GENOME/GRCh37/genome.fa \
    --input GTFGeneModel --input_file $GTF \
    --output CopyNumberAnalysis --output_file $prefix\_FACETS_$data\_copy_number_analysis \
    --dattr FACETSModelSelectionFile.model_selection=adhoc \
    --mattr facets_gbuild=hg19 \
    --mattr cn_header=tcn.em \
    --mattr total_cn_header=tcn.em \
    --mattr minor_cn_header=lcn.em \
    --mattr cn_header2=CNt \
    --mattr total_cn_header2=CNt \
    --mattr minor_cn_header2=B \
    --mattr discard_chrom_with_prefix=MT \
    --run
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

