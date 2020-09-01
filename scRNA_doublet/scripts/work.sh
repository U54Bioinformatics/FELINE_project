
#################run command below to reproduce the results########
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_QC_R3.6/lib/R/library
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/scRNA_QC/lib/python3.6/site-packages/

echo "predict doublet cells with scrublet"
#The input file for this analysis is a seurat obj that have nearly all cells from cell ranger (*_10x_Seurat_2kgenes_vst_cc.raw.rds).
#Output1: FEL011046.scrublet_out.csv contains scrublet prediction of all cells from all patients
#Output2: FEL011046_Raw_individual_results contains scrublet prediction of each individual patients
sbatch --array 1-35 scrublet_QC.sh

echo "merge prediction"
python scrublet_merge_meta.py --list patients.list --output FEL011046
