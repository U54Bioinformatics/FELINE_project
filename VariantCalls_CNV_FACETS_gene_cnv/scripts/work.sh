ln -s ../VariantCalls/FELINE_FACETS_runmerge_copy_number_analysis.txt/seg_files/ ./
cp ../VariantCalls_Maftools/FELINE_FACETS_tumor_purity.txt .
cp ../VariantCalls_phyloWGS/Cancer_genes.CDKi_resistent.txt ./
cp ../Clinical_data//FELINE_patient_1_46.clinical_data.txt ./i

echo "manual selected v1"
cp -R ../VariantCalls/FELINE_FACETS_runmerge_copy_number_analysis.txt/seg_files/ seg_files_20200803_v1

echo "updated difficult ones with new parameters (larger cval), 20200804 v2"
sbatch --array 1-2 FACETS_gene_cnv_reformat_seg.sh
paste ../VariantCalls_Gene_CNV_FACETS_calling/FELINE_FACETS_runmerge.tumor_purity.txt FELINE_FACETS_tumor_purity.txt | cut -f1,2,3,7 > FELINE_FACETS_tumor_purity.20200819.txt
cp FELINE_FACETS_tumor_purity.20200819.txt FELINE_FACETS_tumor_purity.txt
#######################commands below are to reproduce the results ###################################
echo "FACETS segment results"
cp -R ../VariantCalls_Gene_CNV_FACETS_calling/FELINE_FACETS_runmerge_copy_number_analysis/seg_files/ ./

echo "call gene cnv for each sample"
# require four files
# FACETS cnv segment: seg_files 
# Driver gene list:   Cancer_genes.drivers_resistent.list
# Gene coordinate:    hg19.RefSeq.NM_pos_unique_sort.bed
# Sample ploidy:      FELINE_FACETS_tumor_purity.txt
# Sample clinic:      FELINE_patient_1_46.clinical_data.txt 
sbatch --array 1-24 FACETS_gene_cnv_run.sh
mkdir FELINE_FACETS_gene_cnv
mv FEL*.facets.* FELINE_FACETS_gene_cnv/

echo "merge gene cnv from all patients samples"
# generate three tables:
# gene cnv meta table: FELINE_FACETS.gene_cnv_call.short_list.txt
# sample vs. gene CNV call table: FELINE_FACETS.gene_cnv_call.short_list.STable_CNVcall.txt
# sample vs. gene read depth ratio table: FELINE_FACETS.gene_cnv_call.short_list.STable_depthratio.txt
python FACETS_gene_cnv_merge_all_sample.py --input FELINE_FACETS_gene_cnv --output FELINE_FACETS

