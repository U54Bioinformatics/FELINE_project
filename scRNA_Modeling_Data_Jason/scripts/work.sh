
##################################Run commands below to reproduce the results#########################
echo "prepare input files"
cp ../scRNA_meta/output/FEL001046_scRNA.metadata.clinical.txt ./
cp ../FELINE_clinical/output/FELINE_patient_1_46.clinical_data.txt ./
cp ../scRNA_inferCNV_heatmap/output/FEL011046_infercnv_subclones_final_HMM_infer.subclone.anno.txt ./
# cell cycle prediction is not reliable but included anyway
cp ../../Archive/FELINE_merged_011_040_premrna_Archive/Cellcycle/FEL001046_scRNA.cellcycle.txt ./
ln -s ../scRNA_count_zinbwave_merge_10x_ICELL8/output/FEL001046_Cancer_cells_scRNA.zinbwave.normalized.RDS ./
ln -s ../scRNA_ssGSEA_merge_10x_ICELL8/output/FEL001046_Cancer_cells_scRNA.zinbwave.normalized.ssGSEA_scores.RDS ./

echo "Merge data"
# output is FEL011046_data_gene_pathway.RDS, which include patient meta (clinical setting and response) and key genes and pathways.
# gene_list.selected.txt contains genes that Jason and Andrea selected
# gene_list.all_in_RDS.txt contain genes that are in file FEL001046_Cancer_cells_scRNA.zinbwave.normalized.RDS
# in the last version (20200917), I used gene_list.all_in_RDS.txt to include all genes that are avaliable
ln -s gene_list.all_in_RDS.txt gene_list.txt 
sbatch FEL011046_data_gene_pathway.sh
