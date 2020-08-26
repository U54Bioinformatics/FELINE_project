##################Run command below to reproduce results#################

echo "merge clinical information into a single file"
#KI67 data from clinic: KI67_values_at_three_time_points_01-09-2020.txt
#Tumor size (Burden column) from Jason's model: FELINE_clinical_response.txt
#Clinical response classification from Jason's model: FELINE_clinical_response_5groups.txt
#Clinical setting meta data (response not useful in this file): Patient_1_46.clinical_data.metadata.txt
python Classify_ki67_groups.py --ki67 KI67_values_at_three_time_points_01-09-2020.txt --clinical Patient_1_46.clinical_data.metadata.txt --clinical_all FELINE_clinical_response.txt --clinical_5groups FELINE_clinical_response_5groups.txt
# assign two classes response (responder and nonresponder) for final clinical response
# output is "FELINE_patient_1_46.clinical_data.txt"
cat FELINE_patient_1_46.clinical_data.R | R --slave

