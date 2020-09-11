
##################################Run commands below to reproduce the results#########################
echo "convert variant to tsv format"
cp ../VariantCalls_pyclone/output/FELINE_variant_calls_filtered/*_variant_calls_04_purity40_totalRD20_minread5_minVAF0.05.txt ./
awk '{print "python Prepare_somatic_variant_VEP_tsv_FELINE.py --patient "$1}' FELINE_patients.list > run_vep_tsv.sh
bash run_vep_tsv.sh
head -n 1 FEL013_variant_calls_04_purity40_totalRD20_minread5_minVAF0.05.vep.tsv > FELINE.vep.tsv
cat FEL0*.vep.tsv | grep -v "Chrom" >> FELINE.vep.tsv

echo "convert tsv to maf format"
sbatch run_vaf2maf.sh
# make maf format per patient instead of per sample/biopsy
python Somatic_mutation_patient_maf.py --maf FELINE.vep.maf

echo "FACETS gene cnv"
cp ../VariantCalls_CNV_FACETS_gene_cnv/output/FELINE_FACETS.gene_cnv_call4maftools.txt FELINE_FACETS_copy_number_gene.cnv4maftools.txt
# make gene cnv per patient instead of per sample/biopsy
sed 's/_E//g' FELINE_FACETS_copy_number_gene.cnv4maftools.txt | sed 's/_S//g' | sed 's/_M//g' > FELINE_FACETS_copy_number_gene.cnv4maftools.patient.txt

echo "run maftools"
# analysis per patient
# oncoplot output:   FELINE_patient_02_oncoplots.pdf
# oncogenic pathway: FELINE_patient_14_oncogenic_signaling_pathway.pdf 
sbatch run_maftools.patient.sh 
# analysie per sample
# oncoplot output:   FELINE_timepoint_02_oncoplots.pdf
# oncogenic pathway: FELINE_timepoint_14_oncogenic_signaling_pathway.pdf
sbatch run_maftools.sh

#echo "mutated in 2 or more samples"
#cut -f1,11 laml_geneSummary.txt | sort -k2,2rn| awk '$2>=5' > laml_geneSummary.mutated_in_5.txt

