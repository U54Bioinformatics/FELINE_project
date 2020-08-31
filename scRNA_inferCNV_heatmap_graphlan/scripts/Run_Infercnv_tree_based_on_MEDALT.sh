#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=40:00:00
#SBATCH --output=Run_Infercnv_tree_based_on_MEDALT.sh.%A_%a.stdout
#SBATCH -p abild,all
#SBATCH --workdir=./


start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=3
fi

echo "CPU: $CPU"
echo "N: $N"

export R_LIBS=/home/jichen/software/BETSY/install/envs/Infercnv_heatmap_R3.6/lib/R/library:/home/jichen/software/BETSY/install/envs/R_TREE3.6/lib/R/library
Rscript=/home/jichen/software/BETSY/install/envs/Infercnv_heatmap_R3.6/bin/Rscript

patient=`cat patients_arm_response.list | cut -f1 | head -n $N | tail -n 1`
arm=`cat patients_arm_response.list | cut -f2 | head -n $N | tail -n 1`
response=`cat patients_arm_response.list | cut -f3 | head -n $N | tail -n 1`
patient_id=`cat patients_arm_response.list | head -n $N | tail -n 1 | cut -c1-6 | sed 's/FEL0/P/'`
echo "patient IDs"
echo $patient
echo $arm
echo $response
echo $patient_id

# response code
if [ $response == "Non-responder" ]; then
   response="NR"
fi
if [ $response == "Responder" ]; then
   response="R"
fi
echo $response

echo "generate tree for graphlan"
mkdir MEDALT_tree
if [ ! -e MEDALT_tree/$patient\.tree ]; then
   $Rscript Run_Infercnv_tree_based_on_MEDALT.R $patient 
fi

echo "generate annot file for graphlan"
if true; then
   $Rscript Run_graphlan_annot.R $patient
   #echo -e "title\t"$patient_id" (Arm "$arm", "$response")" > $patient\.title.txt
   #cat $patient\.title.txt
   ## draw ESR1, CDK6, FGFR2, ERBB4, RORA
   #cat $patient\.title.txt annot_FELINE.txt $patient\.annot.txt > $patient\.annot1.txt
   cat annot_FELINE.txt $patient\.annot.txt > $patient\.annot1.txt
fi

if false; then
   $Rscript Run_graphlan_annot.fewer_genes.R $patient
   #echo -e "title\t"$patient_id" (Arm "$arm", "$response")" > $patient\.title.txt
   cat $patient\.title.txt
   ## draw ESR1, CDK6, FGFR2, ERBB4
   cat $patient\.title.txt annot_FELINE.fewer_genes.txt $patient\.annot.txt > $patient\.annot1.txt
fi

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

