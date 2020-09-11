
echo "icell8 infercnv"
cp ../scRNA_inferCNV_ICELL8/outputFEL001010_10x_counts.all.output_dir_plot_100cells/infercnv.observations.txt FEL001010.infercnv.observations.txt
# add Gene.ID to first and convert txt to RDS.
sbatch --array 1 ../scripts/Matrix2RDS.sh

