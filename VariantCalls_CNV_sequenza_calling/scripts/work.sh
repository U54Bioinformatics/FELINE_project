########################Run commends below to reproduce results#############
echo "Sequenza CNA calling"
sbatch betsy_WGS_CNA_Sequenza_FELINE_step1.sh

echo "set best model file"
cp FELINE_cnv_sequenza/models.highest_probability.txt ./FELINE_cnv_sequenza_model_selection.txt

echo  "run copy number analysis to generate gene copy numbers and tumor purity, genome instability index"
sbatch betsy_WGS_CNA_Sequenza_FELINE_step3.sh
cp FELINE_cnv_sequenza_copy_number_analysis.txt/copy_number.by_gene.txt FELINE_cnv_sequenza_copy_number_gene.txt
cp FELINE_cnv_sequenza_copy_number_analysis.txt/copy_number.by_segment.txt FELINE_cnv_sequenza_ChromSegmentSig.txt
cp FELINE_cnv_sequenza_copy_number_analysis.txt/tumor_purity.txt FELINE_cnv_sequenza_tumor_purity.txt
cp FELINE_cnv_sequenza_copy_number_analysis.txt/genome_instability_index.txt FELINE_cnv_sequenza_genome_instability_index.txt


