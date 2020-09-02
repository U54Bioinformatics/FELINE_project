
##################Run command below to reproduce results##################

echo "pyclone results folder based on FACETS and sequenze copy number results"
ln -s ../VariantCalls_pyclone_clonevol/output/ VariantCalls_pyclone_clonevol
echo "best model for pyclone with either FACETS or sequenza copy number results"
cp ../VariantCalls_pyclone_clonevol_merge_figure/input/FELINE_pyclone_model.txt ./
echo "variant folder including filtered variants and annotation"
ln -s ../VariantCalls_pyclone/output/FELINE_variant_calls_filtered/ FELINE_variants
echo "driver genes from cosmid"
#https://cancer.sanger.ac.uk/census, use cosmid cancer genes (hg19)
awk -F"," '{print $1}' Census_allMon_Jul27_18_33_19_2020.csv > Cancer_genes.ID.list

echo "class variants: truncal and acquired mutations"
# truncal mutations are defined as these mutations have ccf >= 70% at both timeppint or mutations in a cluster that have ccf >= 70% at both timeppint
# acquired mutations are defined as these mutations have ccf <= 10% at Day 0 and >= 20% at Day 14 or Day 180 
bash FELINE_pyclone_variant_classification.sh

