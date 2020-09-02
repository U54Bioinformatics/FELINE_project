#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=4:00:00
#SBATCH --output=FELINE_pyclone_variant_classification.sh.%A_%a.stdout
#SBATCH -p fast,all
#SBATCH --workdir=./

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


echo "Classify variant for each patient"
for i in {1..24}
do
    patient=`tail -n 24 FELINE_pyclone_model.txt | cut -f1 | head -n $i | tail -n 1`
    cnv=`tail -n 24 FELINE_pyclone_model.txt | cut -f2 | head -n $i | tail -n 1`
    echo "Rank: $i"
    echo "Read: $patient"
    echo "cnv: $cnv"

    if [ ! -e $patient\.variant_info.txt ]; then
        echo "patient: $patient ..."
        python FELINE_pyclone_variant_classification.py --patient $patient --driver Cancer_genes.ID.list --cnv $cnv
    fi
done

echo "Merge table for all patients"
cat *.variant_info.txt | head -n 1 > FELINE_variant_info.txt
cat *.variant_info.txt | grep -v "patient" >> FELINE_variant_info.txt
python ~/software/bin/txt2xlsx.py --input FELINE_variant_info.txt
grep -v "FALSE" FELINE_variant_info.txt > FELINE_variant_info.drivers.txt
python ~/software/bin/txt2xlsx.py --input FELINE_variant_info.drivers.txt


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

