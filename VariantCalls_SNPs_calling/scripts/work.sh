###################run command below to reproduce the results##############
echo "prepare patient individual files for BETSY"
python Prepare_somatic_variant_FELINE.py --patient FELINE_patients.list

echo "run BETSY to call somatic variants"
sbatch --array 1-24 betsy_WGS_variants_FELINE.sh
