
#######################commands below are to reproduce the results ###################################

echo "prepare count and cell meta"
cp ../scRNA_seurat_integration_ICELL8/output/FEL001010_icell8.cell_metadata.UMAPcluster.txt ./
cp ../scRNA_seurat_integration_ICELL8/output/FEL001010_icell8_gene_symbols.raw.counts.txt FEL001010_10x_gene_symbols.filtered.counts.txt
# add Gene.ID to the first row

echo "run infercnv for 10 icell8 patients"
sbatch Infercnv_pipeline_individual_from_individual_filtered_cellline.sh

