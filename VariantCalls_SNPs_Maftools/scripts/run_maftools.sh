#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --output=run_maftools.sh.%A_%a.stdout
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

# prefix for output
prefix=FELINE_timepoint
# gene cnv FACETS, WES
gene_cnv=FELINE_FACETS_copy_number_gene.cnv4maftools.txt
# maf
maf=FELINE.vep.maf
# clinic
clinical=FELINE.clinical.txt

echo "merge two timepoints"
# all analyses
#/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/bin/Rscript Somatic_mutation_MAFtools_FELINE.R FELINE_timepoint FELINE.vep.maf FELINE.clinical.txt $gene_cnv
# oncogenic_signaling_pathway
/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/bin/Rscript Somatic_mutation_MAFtools_FELINE.oncogenic_signaling_pathway.R $prefix $maf $clinical $gene_cnv
# oncoplot
## without cnv
/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/bin/Rscript Somatic_mutation_MAFtools_FELINE.oncoplot.R $prefix\_no_cnv $maf $clinical $gene_cnv
## with cnv
/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/bin/Rscript Somatic_mutation_MAFtools_FELINE.oncoplot_cnv.R $prefix\_cnv $maf $clinical $gene_cnv

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

