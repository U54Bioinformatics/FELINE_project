#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=80G
#SBATCH --time=50:00:00
#SBATCH --output=scRNA_RNA_velocity_step1.sh.%A_%a.stdout
#SBATCH -p all,abild
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh


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

sample=`cat samples.list | head -n $N | tail -n 1`
bam_dir=FEL011046_cell_ranger.mRNA
#GTF=/home/jichen/Projects/Database/Reference/refdata-cellranger-hg19-3.0.0/genes/genes.gtf
GTF=/home/jichen/Projects/Database/Reference/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
TenXfolder=$bam_dir/$sample/
bam_file=$TenXfolder/outs/possorted_genome_bam.bam
out_dir=$TenXfolder/velocyto
loom_file=$TenXfolder/velocyto/$sample\.loom

#obj_seurat=VG_CAMA1_10x_Seurat_2kgenes_vst_cc.rds
#obj_velocity=$sample\_velocity.rds

echo "Step1: Analyze spliced and unspliced reads and genrate a loom file for velocityR"
export PATH=$PATH:/home/jichen/software/BETSY/install/envs/scRNA_velocyto/bin/
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/scRNA_velocyto/lib/python3.6/site-packages/
module load samtools

if [ ! -e $fq_dir/$sample/velocyto/$sample\.loom ]; then
   echo "Run velocyto to generate loom file"
   #Run velocyto tool on cell ranger output folder
   velocyto run10x --samtools-threads $CPU --samtools-memory 500 $TenXfolder $GTF
   #Run velocyto tool on cell ranger (or other) bam files
   #velocyto run --samtools-threads $CPU --samtools-memory 500 -o $out_dir $bam_file $GTF
fi 

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

