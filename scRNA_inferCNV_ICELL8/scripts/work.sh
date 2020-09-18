
#######################commands below are to reproduce the results ###################################

echo "prepare count and cell meta"
# 10x or icell8 does not matter here. Just use 10x to choose a script to run in next step
cp ../scRNA_seurat_integration_ICELL8/output/FEL001010_icell8.cell_metadata.UMAPcluster.txt FEL001010_10x_cell_metadata.UMAPcluster.txt
cp ../scRNA_seurat_integration_ICELL8/output/FEL001010_icell8_gene_symbols.counts.txt FEL001010_10x_gene_symbols.filtered.counts.txt
# remove Gene.ID from the first row if there is one

echo "run infercnv for 10 icell8 patients"
sbatch Infercnv_pipeline_individual_from_individual_filtered_cellline.sh

