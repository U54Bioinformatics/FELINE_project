##########################Run commends below to reproduce results################
### Mutation signature using mutational pattern
echo "Prepare mutation file for 24 patients" 
mkdir FELINE_WES_variant_tab
cp ../VariantCalls_pyclone/FEL*_variant_calls_04_purity40_totalRD20_minread5_minVAF0.05.tab ./FELINE_WES_variant_tab

echo "Run mutational pattern"
sbatch run_Mutationalpattern.sh


