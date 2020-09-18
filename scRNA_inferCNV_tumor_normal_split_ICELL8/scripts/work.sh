#######################commands below are to reproduce the results ###################################
# data from folder "scRNA_inferCNV_tumor_normal_split_10x" is a better example to run tumor-normal split and subclone analysis. 

echo "create meta file for epithelial cells and nonepithelial cells to run tumor-normal split"
# This is how I split tumor-normal for icell8 data after several testing. It should be easier to just run this on all cells. The results should be similar.
head -n 1 FEL001010_cell_metadata.raw.txt > FEL001010_cell_metadata.epi.txt
grep "Epithelial" FEL001010_cell_metadata.raw.txt >> FEL001010_cell_metadata.epi.txt
grep -v "Epithelial" FEL001010_cell_metadata.raw.txt > FEL001010_cell_metadata.noepi.txt

echo "infercnv profile"
# need to remove "Gene.ID" from first row
cp ../scRNA_inferCNV_heatmap/output/FEL001010.infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt ./
cp ../scRNA_inferCNV_heatmap/output/FEL001010.infercnv.observations.txt ./


echo "tumor-normal split for epithelial cells"
# I run this by using epithelial and nonepithelial cells. It does not need to be accurate. Cells were called if has CNA or noCNA. Then the CNA and noCNA call are merged afterward.
ln -s FEL001010_cell_metadata.epi.txt FEL001010_cell_metadata.txt
sbatch infercnv_heatmap_individual_epi.sh
# after runing the scripts, check figures "FEL001010_HMM_infer_k*.png" to select the best split and use "FEL001010_HMM_infer.clust*.anno.txt" to extract the Cluster annotation. In this case, 2 clusters were used, one for tumor (Cluster2) and one for normal. Original output is in folder: FEL001010_epithelial_cells. Rerun results is slightly different and in folder: FEL001010_epithelial_cells_rerun  
python scRNA_Seurat_05individual_from_individual_add_infercnv_tumor_split2cluster.py --infercnv_predict FEL001010_HMM_infer.clust2.anno.txt --meta_cluster FEL001010_cell_metadata.epi.txt --tumor_cluster Cluster2

echo "tumor-normal split for nonepithelial cells"
ln -s FEL001010_cell_metadata.epi.txt FEL001010_cell_metadata.txt
sbatch infercnv_heatmap_individual_noepi.sh
# split tumor and normal based on infercnv heatmap. In this case, 6 clusters were used, 3 for tumor (Cluster1,3,4) and three for normal. Original output is in folder: FEL001010_nonepithelial_cells. Rerun results is slightly different and in folder: FEL001010_nonepithelial_cells_rerun
python scRNA_Seurat_05individual_from_individual_add_infercnv_tumor_split2cluster.py --infercnv_predict FEL001010_HMM_infer.clust6.anno.txt --meta_cluster FEL001010_cell_metadata.noepi.txt --tumor_cluster Cluster1,Cluster3,Cluster4

echo "merge tumor-normal split anno"
cat FEL001010_epithelial_cells/FEL001010_cell_metadata.epi.infercnv_tumor_normal.txt FEL001010_nonepithelial_cells/FEL001010_cell_metadata.noepi.infercnv_tumor_normal.txt | head -n 1 > FEL001010_cell_metadata.infercnv_tumor_normal.txt
cat FEL001010_epithelial_cells/FEL001010_cell_metadata.epi.infercnv_tumor_normal.txt FEL001010_nonepithelial_cells/FEL001010_cell_metadata.noepi.infercnv_tumor_normal.txt | grep -v "Cell.ID" >> FEL001010_cell_metadata.infercnv_tumor_normal.txt

echo "subclone"
# I did not run subclone anlaysis for icell8 as too few patients have cancer cells avaliable pre- and post treatment 

