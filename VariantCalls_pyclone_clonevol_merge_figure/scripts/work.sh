echo "data and script"
cp ../VariantCalls_pyclone_clonevol/FEL022_pyclone_analysis_FACETS_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.RData ./
cp ../VariantCalls_pyclone_clonevol/FELINE_clonevol_ccf.list ./
cp ../VariantCalls_pyclone_clonevol/FELINE_Step3_run_clonevol.sh ./FELINE_Step4_run_fishplot.sh
cp ../Clinical_data/FELINE_patient_1_46.clinical_data.txt ./

#######################run commend below to reproduce results
echo "cp fishplot data for best pyclone model/cnv: FEL021 and FEL029 failed to generate pyclone results"
cat FELINE_pyclone_model.txt | grep "^FEL" | awk -F"\t" '{print "cp ../VariantCalls_pyclone_clonevol/"$1"_pyclone_analysis_"$2"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.RData ./"}' > FELINE_pyclone_model.cp.sh
bash FELINE_pyclone_model.cp.sh

echo "tumor size changes"
bash FELINE_clinical_response.sh

echo "run fishplot, need to finish FEL025"
ls FEL0*_pyclone_analysis_*20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.RData | sed 's/clonevol_figure4_fishplot.RData/txt/' > FELINE_clonevol_ccf.list
bash FELINE_Step4_run_fishplot.sh

echo "merge fishplot into one figure"
cut -f1,2,3,4,5 FELINE_patient_info.sorted.txt > FELINE_patient_info.sorted.revised.txt
bash FELINE_fishplot.sh

