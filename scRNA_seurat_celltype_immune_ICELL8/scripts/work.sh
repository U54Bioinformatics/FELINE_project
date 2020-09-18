
#######################commands below are to reproduce the results ###################################

echo "extrat candidate immune cells (meta and count), this step can be skipped for icell8" 
cp ../scRNA_seurat_integration_ICELL8/output/FEL001010_icell8_gene_symbols.counts.txt ./
cp ../scRNA_seurat_integration_ICELL8/output/FEL001010_icell8.cell_metadata.UMAPcluster.txt ./
grep "Immune" ../scRNA_seurat_celltype_annotation_ICELL8/output/FEL001010_icell8.cell_metadata.UMAPcluster.revised_anno.txt | cut -f1 > FEL001010_icell8_anno.Immune.txt
python Extract_cell_meta.py --cell FEL001010_icell8_anno.Immune.txt --meta FEL001010_icell8.cell_metadata.UMAPcluster.txt 
python Extract_cell_count.py --cell FEL001010_icell8_anno.Immune.txt --count FEL001010_icell8_gene_symbols.counts.txt

echo "cluster immune cells, this step can be skipped for icell8"
ln -s FEL001010_icell8_gene_symbols.counts.txt.subset.txt FEL001010_icell8.immune.txt
ln -s FEL001010_icell8.cell_metadata.UMAPcluster.txt.subset.txt FEL001010_icell8.immune.metadata.txt
sbatch scRNA_Seurat_count2umap.sh

echo "Finalize immune cell annotation based on singleR annotation"
# try to use seurat cluster from above analysis. but finally choose to use singleR annotation alone.
# input 1: FEL001010_icell8.immune.metadata.txt
# output 1: FEL001010_icell8.immune.metadata.revised_anno.txt
# in output file two new columns were added (Encode_main_cluster     Encode_main_cluster_revised) to define immune cell classfication. Encode_main_cluster is broader (T cells) and Encode_main_cluster_revised is more detailed (CD4+ or CD8+ T cells). 
python scRNA_add_cell_type_immune_icell8.py --meta FEL001010_icell8.immune.metadata.txt

