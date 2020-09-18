#######################commands below are to reproduce the results ###################################

echo "I use FEL046P117 as example to show how to run this analysis. The final results for all patients are in scRNA_inferCNV_heatmap"
cp ../scRNA_inferCNV_10x/output/FEL011046_infercnv_result_folder/FEL046P117_10x_counts.all.output_dir_plot_100cells/FEL046P117_* ./

echo "run infercnv on all cells. the results can be used to split tumor and normal cells"
# output 1: FEL046P117_HMM_infer_k*.png
# output 2: FEL046P117_HMM_infer.clust*.anno.txt
# in this case, 9 clusters were used and Cluster1,Cluster2,Cluster3,Cluster4,Cluster5,Cluster6,Cluster9 were found be tumor clusters
sbatch --array 1 infercnv_heatmap_individual.sh
# results for all 35 10x patients are in "FEL011046_infercnv_subclones"

echo "split tumor and normal"
# output: FEL046P117_cell_metadata.infercnv_tumor_normal.txt;
# column "Infercnv_CNA" indicate whather a cell has CNA or noCNA
# column "Infercnv_split" indicate whether tumor cell or normal cell  
python scRNA_Seurat_05individual_from_individual_add_infercnv_tumor_split2cluster.py --infercnv_predict FEL046P117_HMM_infer.clust9.anno.txt --meta_cluster FEL046P117_cell_metadata.txt --tumor_cluster Cluster1,Cluster2,Cluster3,Cluster4,Cluster5,Cluster6,Cluster9
# run_tumor_normal_split.sh records tumor clusters for each patient 
# results for all 35 10x patients are in "FEL011046_revised_anno"

echo "extract meta file for subclone analysis"
# use Infercnv_split to extract cells that are tumor
# output: FEL046P117_cell_metadata_subclone.txt
python infercnv_heatmap_get_cell_metadata.py --meta FEL046P117_cell_metadata.infercnv_tumor_normal.txt --subclone
# results for all 35 10x patients are in "FEL011046_cell_metadata_subclone"

echo "run subclone analysis"
# output 1: FEL046P117_HMM_infer.subclone_k*.png
# output 2: FEL046P117_HMM_infer.subclone*.anno.txt 
# run 5, 6, 7 clusters
# 5 clusters are chosen.
sbatch --array 5-7 infercnv_heatmap_individual_subclone46.sh
# results for all 35 10x patients are in "FEL011046_infercnv_subclones" or "FEL011046_infercnv_subclones_final". The later has only the chosen one for each patient
