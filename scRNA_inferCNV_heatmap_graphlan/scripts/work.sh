
#################run command below to reproduce the results########
echo "prepare input files"
# data file prepared for Jason including patient meta information and gene expression
cp ../../scRNA_Modeling_Data_Jason/output/FEL011046_data_gene_pathway.RDS
# patient id/arm/response
# FEL011P138	B	Non-responder
# FEL012P101	B	Non-responder
cp ../FELINE_clinical/output/FELINE_patient_1_46.clinical_data.txt ./
cut -f1,4,11 FELINE_patient_1_46.clinical_data.txt | tail -n 36 | grep -v "FEL018P110" > patients_arm_response.list
# infercnv copy number profile reformated with MEDALT
mkdir MEDALT_cnv
cp ../scRNA_inferCNV_heatmap_MEDALT/output/FEL011046_infercnv_MEDALT_results/FEL0*.HMMi6.final.CNV.txt MEDALT_cnv/
cp ../scRNA_inferCNV_heatmap_MEDALT/output/FEL011046_infercnv_MEDALT_results/FEL0*.HMMi6.final.CNV.names.txt MEDALT_cnv/
# template config file for graphlan
annot_FELINE.txt
annot_FELINE.fewer_genes.txt

echo "prepare config file and phylogenetic tree based on infercnv copy number profile for 35 patients"
# input are "patients_arm_response.list", "MEDALT_cnv" and "FEL011046_data_gene_pathway.RDS"
# outout 1: a folder "MEDALT_tree" including trees for 35 patients
# output 2: config files "FEL0*annot1.txt" for graphlan to generate circos plot 
sbatch --array 1-35 Run_Infercnv_tree_based_on_MEDALT.sh

echo "generate circos plot with graphlan"
# input are "patients_arm_response.list", "FEL0*annot1.txt", "MEDALT_tree/*.tree", 
# output are "FEL0*.tree.pdf/png"
sbatch --array 1-35 Run_graphlan_tree.sh



