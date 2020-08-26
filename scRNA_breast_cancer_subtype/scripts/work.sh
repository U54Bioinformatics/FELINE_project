
echo "prepare CPM for each patient"
ln -s ../scRNA_count_CPM/output/FEL011046_10x_FEL0*.CPM.txt ./
ls *.CPM.txt | sed 's/_gene_symbols.CPM.txt//g' > FEL011046.patients.list

echo "prepare cell type annotation and a gene list that we will use to summary genefu results"
# output is "FEL011046.cell_type.txt"
#Cell.ID Sample  Celltype1       Celltype2
#FEL011_M_AAACGAACACAAGTGG       FEL011P138_M    Cancer cells    Cancer cells
#FEL011_M_AAACGCTGTTAAGAAC       FEL011P138_M    Cancer cells    Cancer cells
#FEL011_M_AAAGGTATCAGCCTCT       FEL011P138_M    Macrophages     Macrophages
# output is "FEL011046.gene_id.txt"
#Gene.ID_1       Gene.ID_2       Gene.Symbol
#ENSG00000243485 107985730       MIR1302-2HG
#ENSG00000237613 645520  FAM138A
#ENSG00000186092 79501   OR4F5
sbatch Run_celltype_file.sh

echo "run genefu for each patient"
# many output files for each patient
# these two files summarize number/percent of each intrinsic subtype in each patient
# 1. FEL011046_10x_FEL046_gene_symbols.genefu.sample_sum.txt
# 2. FEL011046_10x_FEL046_gene_symbols.genefu.sample_pc.txt
sbatch --array 1-35 Run_genefu_scRNA.sh

echo "merge patients and plot heatmap for subtype across timepoint"
# merge "*.genefu.sample_sum.txt" and "*.genefu.sample_pc.txt" files
# output are FEL011046.patients_gene_symbols.genefu.sample_sum.txt and FEL011046.patients_gene_symbols.genefu.sample_pc.txt
python scRNA_Seurat_05individual_from_individual_merge_matrix.py --list FEL011046.patients.list
# prepare matrix for heatmap 
# output is FEL011046.subtype_heatmap_matrix.data.txt
python Prepare_subtype_per_patient_heatmap.py --subtype FEL011046.patients_gene_symbols.genefu.sample_pc.txt --output FEL011046.subtype_heatmap_matrix.data.txt
# plot heatmap
# output is FEL011046.subtype_heatmap_matrix.heatmap.pdf/png 
cp ../FELINE_clinical/output/FELINE_patient_1_46.clinical_data.txt ./
bash FEL011046.subtype_heatmap_matrix.sh


echo "compare change across timepoint"
python Prepare_subtype_per_sample.py --patient_data FEL011046_patients_gene_symbols.genefu_filter.sample_sum.txt --patient_responder ../Clinical_data/FELINE_patient_1_46.clinical_data.txt --subtype Basal --output FEL011046


#the published SSPs (SSP203, SSP2006, and PAM50 composed of 500, 306, and 50 genes, respectively) and SCMs (SCMOD1, SCMOD2, and SCMGENE, composed of 726, 663, and 3 genes, respectively).
#scmod1, 726 genes 

