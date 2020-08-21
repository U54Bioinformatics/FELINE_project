echo "data and script"
cp ../VariantCalls/betsy_WGS_CNA_FACETS_FELINE_step1.sh .
cp ../VariantCalls/betsy_WGS_CNA_FACETS_FELINE_step2.sh .
cp ../VariantCalls/betsy_WGS_CNA_FACETS_FELINE_step3.sh .
ln -s ../VariantCalls/FELINE_bam ./
cp -R ../VariantCalls/FELINE_FACETSResults_runmerge/ ./
cp -R ../VariantCalls/FELINEwes_FACETSResults_20_25_36_v2/facets.cval_1000.ndepth_35.nbhd_250.nhet_15/ ./FELINE_FACETSResults_runmerge/
cp -R ../VariantCalls/FELINEwes_FACETSResults_ontarget_v1/facets.cval_* FELINE_FACETSResults_runmerge/

##################################Run commands below to reproduce the results#########################
echo "run FACETS with multiple parameters"
sbatch betsy_WGS_CNA_FACETS_FELINE_step1.sh

echo "skip betsy_WGS_CNA_FACETS_FELINE_step2.sh, which do similar things as step3"
#what I did is to:
#1. move FACETS result folder (for example facets.cval_1000.ndepth_50.nbhd_250.nhet_15) into "FELINE_FACETSResults_runmerge"
#2. manually select best model and add the parameters into "FELINE_FACETSResults_runmerge/FELINE_FACETS_runmerge_model_selection.txt"

echo "run BETSY to summarize gene level based on FACETS best model"
cp FELINE_FACETSResults_runmerge/FELINE_FACETS_runmerge_model_selection.txt ./
sbatch betsy_WGS_CNA_FACETS_FELINE_step3.sh
cp FELINE_FACETS_runmerge_copy_number_analysis/copy_number.by_gene.txt FELINE_FACETS_runmerge.copy_number.by_gene.txt
cp FELINE_FACETS_runmerge_copy_number_analysis/copy_number.by_segment.txt FELINE_FACETS_runmerge.copy_number.by_segment.txt
cp FELINE_FACETS_runmerge_copy_number_analysis/tumor_purity.txt FELINE_FACETS_runmerge.tumor_purity.txt
