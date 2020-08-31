
##################Run command below to reproduce results##################

echo "input files: infercnv profile and subclone classification"
cp -R ../scRNA_inferCNV_heatmap/output/FEL011046_infercnv_subclones_final ./
cp -R ../scRNA_inferCNV_heatmap/output/FEL011046_infercnv_HMMi6.observations ./

echo "input file: data file prepared for Jason including patient meta information and gene expression"
cp ../scRNA_Modeling_Data_Jason/output/FEL011046_data_gene_pathway.RDS ./

echo "generate scRNA graph based tree using MEDATL and plot with igraph"
# run all cells is very slow and has memory issues
# use random 500 cells instead
# need to change false to true in "Run_MEDALT.sh" to control the steps to run
sbatch --array 1-35 Run_MEDALT.sh


