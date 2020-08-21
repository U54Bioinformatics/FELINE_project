#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=4:00:00
#SBATCH --output=FACETS_gene_cnv_run.sh.%A_%a.stdout
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

if [ $N == 1 ] & [ $1 -gt 1 ];then
   N=$1
   echo "N: $N"
fi

#sample=FEL031
sample=`cat FELINE_patients.list | head -n $N | tail -n 1`

if [ $sample == 'FEL044' ] || [ $sample == 'FEL045' ] || [ $sample == 'FEL046' ] || [ $sample == 'FEL036' ]; then
    if [ ! -e $sample\_M.facets.gene_cnv_call.short_list.txt ];then
        python FACETS_gene_cnv.py --facets seg_files/$sample\_S.seg --sample $sample\_S
        python FACETS_gene_cnv.py --facets seg_files/$sample\_M.seg --sample $sample\_M
    fi
else
    if [ ! -e $sample\_E.facets.gene_cnv_call.short_list.txt ];then
        python FACETS_gene_cnv.py --facets seg_files/$sample\_S.seg --sample $sample\_S
        python FACETS_gene_cnv.py --facets seg_files/$sample\_E.seg --sample $sample\_E
    fi
fi


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

