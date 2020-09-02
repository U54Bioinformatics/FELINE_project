##################Run command below to reproduce results##################
### Run 24 patients with sequenza and facets cnv

export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol_R/lib/R/library

echo "cp sequenza/FACETS results folder"
bash FELINE_Step1_cp_results_folder.sh

echo "prepare infile for clonevol"
ls -d FEL0*_pyclone_analysis_* | grep -v "ccf" | grep -v "vaf" > FELINE_clonevol_prepare_table.list
bash FELINE_Step2_prepare_clonevol_infile.sh

echo "run for all patients"
ls *.clonevol_ccf.txt > FELINE_clonevol_ccf.list
bash FELINE_Step3_run_clonevol.sh
mv *.pdf FELINE_clonevol_sequenza_FACETS_ccf

