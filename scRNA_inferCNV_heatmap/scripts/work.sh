#######################commands below are to reproduce the results ###################################


echo "icell8 infercnv observation for Supplementary Figure 1 infercnv heatmap"
# output 1: FEL001010.infercnv.observations.txt, matrix of infercnv observation from 10 icell8  patients
# output 2: FEL001010.infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt, matrix of infercnv HMM observation from 10 icell8  patients
cp ../scRNA_inferCNV_ICELL8/output/FEL001010_10x_counts.all.output_dir_plot_100cells/infercnv.observations.txt FEL001010.infercnv.observations.txt
cp ../scRNA_inferCNV_ICELL8/output/FEL001010_10x_counts.all.output_dir_plot_100cells/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt FEL001010.infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt
# add Gene.ID as a header to first row in FEL001010.infercnv.observations.txt
# convert txt to RDS
sbatch --array 1 Matrix2RDS.sh

echo "10x infercnv observation for Figure 1 infercnv heatmap"
# output is FEL011046.infercnv.observations.RDS, a merged matrix of infercnv observation from 35 10x patients 
mkdir FEL011046_infercnv.observations
cp ../scRNA_inferCNV_10x/output/FEL011046_infercnv_result_folder/FEL0*/FEL0*_infercnv.observations.txt FEL011046_infercnv.observations/
ls FEL011046_infercnv.observations/FEL0*.txt | sed 's/infercnv\.observations\.//' > file_infercnv.txt
sbatch Normalize_count_by_zinbwave_merge_output.sh
# convert txt to RDS
sbatch --array 1 Matrix2RDS.sh

echo "10x infercnv HMM observation for scRNA phylogeny (MEDALT and graphlan)"
mkdir FEL011046_infercnv_HMMi6.observations
cp ../scRNA_inferCNV_10x/output/FEL011046_infercnv_result_folder/FEL0*/*_infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt FEL011046_infercnv_HMMi6.observations/
# no need to merge HMM observation as we use as individual patient

echo "10x infercnv final subclone"
cp -R ../scRNA_inferCNV_tumor_normal_split_10x/output/FEL011046_infercnv_subclones_final ./
ls FEL011046_infercnv_subclones_final/FEL0* | sed 's/_HMM_infer.subclone.anno.txt//' > FEL011046_infercnv_subclones_final.patient.list
python infercnv_merge_subclone.py --list FEL011046_infercnv_subclones_final.patient.list --output FEL011046_infercnv_subclones_final

echo "10x infercnv tumor-normal split"
cp -R ../scRNA_inferCNV_tumor_normal_split_10x/output/FEL011046_revised_anno/ ./
ls FEL011046_revised_anno/FEL0* | sed 's/_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.infercnv_tumor_normal.txt//' | grep -v ".txt" > FEL011046_revised_anno.patient.list
python infercnv_merge_meta.py --list FEL011046_revised_anno.patient.list --output FEL011046
mv FEL011046_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.infercnv_tumor_normal.txt FEL011046_cell_metadata.infercnv_CNA.txt


