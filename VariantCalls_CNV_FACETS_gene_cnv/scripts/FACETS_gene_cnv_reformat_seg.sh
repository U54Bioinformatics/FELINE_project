#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=4:00:00
#SBATCH --output=FACETS_gene_cnv_reformat_seg.sh.%A_%a.stdout
#SBATCH -p abild,all
#SBATCH --workdir=./

export PYTHONPATH=/home/jichen/software/BETSY/install/envs/phyloWGS/lib/python2.7/site-packages/
export LD_LIBRARY_PATH=/home/jichen/software/BETSY/install/envs/phyloWGS/lib/:$LD_LIBRARY_PATH

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
echo "Command line N: $1"

if [ $N == 1 ] & [ $1 > 1 ];then
   N=$1
   echo "N: $N"
fi

#sample=FEL036
sample=`cat FELINE_FACETS_parameters.list | grep -v "^#" | head -n $N | cut -f1 | tail -n 1`
parameter=`cat FELINE_FACETS_parameters.list | grep -v "^#" | head -n $N | cut -f2 | tail -n 1`
echo "$sample\t$parameter"

source_dir="../VariantCalls/FELINEwes_FACETSResults_ontarget_v1"
target_dir="./seg_files"
infile=$source_dir/$parameter/$sample\.seg
outfile=$target_dir/$sample\.seg
Rscript FACETS_gene_cnv_reformat_seg.R $infile $outfile
    

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

