
#######################commands below are to reproduce the results ###################################
# Note 1. infercnv results were generate initially when we have patient data arrived. This is a computational intensive process. I did not rerun these analysis. Here I used two patients to show how I ran infercnv for patient with fewer cells (<5000: FEL042P401 has 3k) and patients with large number of cells (>5000: FEL043P202 has 13k). The whole analysis can be reproduced using these two approaches.
# Note 2. I ran infercnv with filtered cells. But I think this should be done with all cells then filtered afterward. Because we don't want to rerun this again due to change of filtering in the analysis. Low quality may affect subclone classification but not infercnv run. If run with all cells a filter is needed before subclone classification.
echo "Input files"
# input 1: hg19.RefSeq.NM_pos_unique_sort.txt, hg19 gene coordinate
# input 2: FELINE.Normal.annot.txt, annotation of 500 random normal cells from FELINE used as reference
# input 3: FELINE.Normal.counts.txt, count of 500 random normal cells from FELINE used as reference
# input 4: FEL041P137_10x_cell_metadata.UMAPcluster.txt, cell meta data only use col1 (cell id) and 5 (sample id).
# input 5: FEL041P137_10x_gene_symbols.filtered.counts.txt, count data (first column is geneid but without header)
## generate meta and count files
ln -s ../scRNA_seurat_integration_10x/input/FEL0*_10x_Seurat_2kgenes_vst_cc.raw.rds ./
ls FEL*.raw.rds| sed 's/_10x_Seurat_2kgenes_vst_cc.raw.rds//' > patients.list
sbatch --array 1-35 scRNA_Seurat_export_count_matrix.sh

echo "Run FEL042P401 with 3k cells"
# create a patients.infercnv.list file with patient id for SLURM array job
# patients.infercnv.list
# FEL042P401
# run only 1 array job as we have fewer cells. we want infercnv to run all cells in one run.
sbatch --array 1 Infercnv_pipeline_individual_from_individual_filtered_cellline.sh
echo "Run FEL043P202 with 13k cells"
# split cells into three chunks each has 4.4k cells; 
python Infercnv_pipeline_individual_from_individual_split_infercnv_input.py --anno FEL043P202_10x_cell_metadata.UMAPcluster.txt --count FEL043P202_10x_gene_symbols.filtered.counts.txt
# create a patients.infercnv.list file with chunk id for SLURM array job
# patients.infercnv.list
# FEL043P202a
# FEL043P202b
# FEL043P202c
# run three 3 array job each for a infercnv analysis on one chunk of cells
sbatch --array 1-3 Infercnv_pipeline_individual_from_individual_filtered_cellline.sh
# merge infercnv results after runing on splitted files
# --cell 500 is the parameter used in line 86 of Infercnv_pipeline_individual_from_individual_filtered_cellline.sh 
# min_cells_per_gene, control resolution of infercnv; use 100 or 500 and choose best results based on heatmap
python Infercnv_pipeline_individual_from_individual_merge_infercnv_obs.py --subsample FEL043P202a,FEL043P202b,FEL043P202c --cell 500


echo "Collect raw infercnv results"
# I copied the data from the original place into a master folder: FEL011046_infercnv_result_folder.
# FEL011 to FEL024
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL011P138_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL012P101_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL013P102_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL014P103_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL015P105_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL016P107_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL017P109_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL019P113_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL020P115_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_500cells/FEL021P118_10x_counts.all.output_dir_plot_500cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL022P124_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL023P121_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
cp -R /net/isi-dcnl/ifs/user_data/abild/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_024_premrna/Infercnv/FEL011024_infercnv_filtered_control_FELINEnormal_all_run_100cells/FEL024P123_10x_counts.all.output_dir_plot_100cells/ ./FEL011046_infercnv_result_folder/
# FEL025 to FEL028
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_025_028_premrna/InferCNV/FEL025P142_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_025_028_premrna/InferCNV/FEL026P140_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_025_028_premrna/InferCNV/FEL027P125_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_025_028_premrna/InferCNV/FEL028P131_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
# FEL029 to FEL040
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL029P134_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL030P135_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL031P143_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL032P145_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL033P129_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL034P601_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL035P119_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL036P132_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL037P703_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL038P403_10x_counts.all.output_dir_plot_100cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL039P201_10x_counts.all.output_dir_plot_100cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_029_040_premrna/InferCNV/FEL040P701_10x_counts.all.output_dir_plot_500cells ./FEL011046_infercnv_result_folder/
# FEL041 to FEL043
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_041_043_premrna/InferCNV/FEL041P137_10x_counts.all.output_dir_plot_100cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_041_043_premrna/InferCNV/FEL042P401_10x_counts.all.output_dir_plot_100cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_041_043_premrna/InferCNV/FEL043P202_10x_counts.all.output_dir_plot_100cells ./FEL011046_infercnv_result_folder/
# FEL044 to FEL046
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_044_046_premrna/InferCNV/FEL044P116_10x_counts.all.output_dir_plot_100cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_044_046_premrna/InferCNV/FEL045P122_10x_counts.all.output_dir_plot_100cells ./FEL011046_infercnv_result_folder/
cp -R /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_044_046_premrna/InferCNV/FEL046P117_10x_counts.all.output_dir_plot_100cells ./FEL011046_infercnv_result_folder/

echo "rename observation files to FEL045P122_observation.txt format"
ls -d FEL011046_infercnv_result_folder/FEL0* | awk '{print "python Infercnv_rename_files_FEL011046.py --infercnv_outdir "$1}' > Infercnv_rename_files_FEL011046.sh
bash Infercnv_rename_files_FEL011046.sh

## keep only observation txt files
#rm FEL0*/infercnv.*
#rm FEL0*/*_obj
#rm FEL0*/*.dat
#rm FEL0*/*.cell_groupings
#rm -R FEL0*/BayesNetOutput.*
