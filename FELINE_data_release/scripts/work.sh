
```
1.  A table with columns for:
- patient ID  (please use the FELxxx designation rather than the original FELINE IDs, so it's harder for people to match them)
- time point
- any other attributes that you will use in the paper (treatment, ki67 status, outcomes, etc)

2.  The BAM files with the whole exome sequencing data for each patient
The names of the files should be:
<patient_id>_<time_point>.bam

3.  One (or more, if it's too big) tab-delimited text files with the single-cell RNA-Seq.  The names of the cells should be <patient_ID>_<time point>_<barcode>.
```

echo "1. patient meta"
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
cp ../FELINE_clinical/output/FELINE_patient_1_46.clinical_data.txt ./
/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript FELINE_patient_1_46.clinical_data.release.R

echo "2. WES bam files"
#ln -s /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_bam/output/FELINE_bam/ output/
sed 's/FEL0//' FELINE_patients.list | awk '{print "ln -s /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_bam/output/FELINE_bam/FEL0"$1"_G.bam ./output/FELINE_bam/P"$1"_G.bam"}' > link_bam_files.sh
sed 's/FEL0//' FELINE_patients.list | awk '{print "ln -s /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_bam/output/FELINE_bam/FEL0"$1"_S.bam ./output/FELINE_bam/P"$1"_S.bam"}' >> link_bam_files.sh
sed 's/FEL0//' FELINE_patients.list | awk '{print "ln -s /home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_bam/output/FELINE_bam/FEL0"$1"_E.bam ./output/FELINE_bam/P"$1"_E.bam"}' >> link_bam_files.sh
# edit P36_E to P36_M, P44_E to P44_M, P45_E to P45_M, P46_E to P46_M
bash link_bam_files.sh

echo "3. scRNA count files"
# icell8
cut -f1 ../Archive/Normalize_count_icell8/FEL001010_Cancer_cells_icell8.metadata.txt |grep -v "Cell.ID" > FEL001010_Cancer_cells_icell8.Cell_ID.txt
cp ../Archive/Normalize_count_icell8/FEL001010_Cancer_cells_icell8.count.txt ./
## fix icell8 cell ids
cp ../scRNA_meta/output/FEL001046_scRNA.metadata.clinical.txt ./
python fix_cell_id_icell8.py --input FEL001010_Cancer_cells_icell8.count.txt --meta FEL001046_scRNA.metadata.clinical.txt
sed 's/FEL0/P/g' FEL001010_Cancer_cells_icell8.count.id_fixed.txt > FEL001010_Cancer_cells_icell8.count.release.txt
# 10x
cp ../Archive/Normalize_count_10x/FEL011046_Cancer_cells_10x.Cell_ID.txt ./
cp ../Archive/Normalize_count_10x/FEL011046_Cancer_cells_10x.count.txt ./
## check 10x cell id: good
python check_cell_id_10x.py --input FEL011046_Cancer_cells_10x.Cell_ID.txt --platform 10x
sed 's/FEL0/P/g' FEL011046_Cancer_cells_10x.count.txt > FEL011046_Cancer_cells_10x.count.release.txt


