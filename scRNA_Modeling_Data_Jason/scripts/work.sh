
##################################Run commands below to reproduce the results#########################
echo "prepare input files"
cp ../scRNA_meta/output/FEL001046_scRNA.metadata.clinical.txt ./
cp ../FELINE_clinical/output/FELINE_patient_1_46.clinical_data.txt ./
cp ../scRNA_inferCNV/output/FEL011046_infercnv_subclones_final_HMM_infer.subclone.anno.txt  FEL011046_infercnv_subclones_final_HMM_infer.subclone.anno.txt
# cell cycle prediction is not reliable but included anyway
cp ../Archive/Cellcycle/FEL001046_scRNA.cellcycle.txt ./
ln -s ../scRNA_count_zinbwave_merge_10x_ICELL8/output/FEL001046_Cancer_cells_scRNA.zinbwave.normalized.RDS ./
ln -s ../scRNA_ssGSEA_merge_10x_ICELL8/output/FEL001046_Cancer_cells_scRNA.zinbwave.normalized.ssGSEA_scores.RDS ./

echo "Merge data"
# output is FEL011046_data_gene_pathway.RDS, which include patient meta (clinical setting and response) and key genes and pathways. 
sbatch FEL011046_data_gene_pathway.sh
