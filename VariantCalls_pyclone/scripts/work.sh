#################run commend below to reproduce the results########
echo "Prepare input file for BETSY to run pyclone"
# bam files, normal_cancer.txt, sample.txt and variant_calls for each patient from "VariantCalls_SNPs_calling" 
ln -s /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_SNPs_calling/output/FEL0* ./
# cnv results from FACETS and sequenza
ln -s /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_CNV_sequenza_calling/output/FELINE_cnv_sequenza ./
ln -s /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_CNV_FACETS_calling/output/FELINE_FACETSResults_runmerge FELINE_FACETSResults
# create FACETS and sequenza folder for each patient
python Prepare_somatic_cnv.py --patient FELINE_patients.list

echo "run mutation filter and pyclone with BETSY for 24 patients"
sbatch --array 1-24 betsy_WGS_pyclone.sh

echo "sum mutation number in pyclone analysis"
bash run_pyclone_mutation_num.sh

echo "estimate tumor burden"
sbatch --array 1-24 betsy_WGS_tumor_burden_FELINE.sh
cat FEL0*_tumor_burden.txt/tumor_mutation_burden.txt | head -n 1 > FELINE_tumor_mutation_burden.txt
cat FEL0*_tumor_burden.txt/tumor_mutation_burden.txt | grep -v "Sample" >> FELINE_tumor_mutation_burden.txt
mv FEL0*_tumor_burden.txt FELINE_tumor_mutation_burden

echo "count INDELs"
python Prepare_somatic_variant_count_INDEL.py --input FELINE_patients.list --variant_folder FELINE_variant_calls_filtered/ > FELINE_WES_INDEL_count.txt

echo "merge tumor purity, variants and tumor burden"
cp ../VariantCalls_CNV_FACETS_gene_cnv/input/FELINE_FACETS_tumor_purity.txt ./

