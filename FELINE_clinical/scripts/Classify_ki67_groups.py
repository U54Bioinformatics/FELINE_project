#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import re
import os
import argparse

def usage():
    test="name"
    message='''


--input: Ki67 file
ID from project managers/VELOS  Correct ID to match with randomization list/treatment arm       Arm     Baseline KI67 (%)       Baseline KI67 original entry    Day 14 KI67 (%) D
001-101 2972-001-101    B       15      0.15    2.24    2.24    5.8     5.8     Letrozole + 600mg Ribociclib (L+Ri)
001-102 2972-001-102    C       12      0.12    0.49    0.49    3.1     3.1     Letrozole + 400mg Ribociclib (L+Rc)
001-103 2972-001-103    A1      5       0.05    0       <1      0.9     0.9     Letrozole + placebo (L+P)


    '''
    print message


def classify_ki67_group(infile, patient_response, patient_response_v3):
    data_ki67 = defaultdict(lambda : str())
    data_v3   = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=0)
    ofile = open("FELINE_all_patients.ki67_groups.txt", 'w')
    print >> ofile, 'Patient.Study.ID\tArm\tKi67_Response\tResponse_v2\tReponse_v3\tReponse_v3b'
    for i in range(file_df.shape[0]):
        trail_id = file_df['ID from project managers/VELOS'][i]
        ki67_d14  = file_df['Day 14 KI67 (%)'][i]
        ki67_d180 = file_df['Surgery KI67 (%)'][i]
        #print file_df.iloc[i,]
        ki67_group = ''
        if float(ki67_d14) >= 0 and float(ki67_d180) >= 0:
            #print "okay"
            #comments from Dr. Qmar and Andrea
            #Sensitive = D14 < 2.7 and Surgery < 2.7
            #Innate resistance: D 14 Ki 67 > 10% ( if we need to expand this group, we can include patients who are >2.7 at D14)
            #Acquired resistance (on therapy). This is the  most interesting group. < 2.7 at D14 but >2.7 sat surgery. 
            if float(ki67_d14) <= 2.7 and float(ki67_d180) <= 2.7:
                ki67_group = 'Sensitive'
            elif float(ki67_d14) <= 2.7 and float(ki67_d180) > 2.7:
                ki67_group = 'Acquired resistance'
            elif float(ki67_d14) > 2.7 and float(ki67_d14) <= 10:
                ki67_group = 'Innate resistance'
            elif float(ki67_d14) > 10:
                ki67_group = 'Innate resistance'
            else:
                ki67_group = 'Unclassified'
        else:
            #print "nodata"
            ki67_group = 'NoData'
        arm = file_df['Arm'][i]
        if arm in ['A1', 'A2']:
            arm = 'A'
        #response = ''
        #responsev3 = ''
        #if patient_response.has_key(trail_id):
        #    response  = patient_response[trail_id][0]
        #else:
        #    response = 'NoData'
              
        response   = patient_response[trail_id][0] if patient_response.has_key(trail_id) else 'NoData'
        responsev3 = patient_response_v3[trail_id][0] if patient_response_v3.has_key(trail_id) else 'NoData'
        responsev3b = patient_response_v3[trail_id][1] if patient_response_v3.has_key(trail_id) else 'NoData'

        print >> ofile, '%s\t%s\t%s\t%s\t%s\t%s' %(trail_id, arm, ki67_group, response, responsev3, responsev3b)
        data_ki67[trail_id] = ki67_group
        data_v3[trail_id] = [responsev3, responsev3b]
    return data_ki67, data_v3

#../../Clinical_Data/FEL011043/FELINE_clinical_response.txt
#Patient.Study.ID        Day     Clinical assessment     axillary_lymph  Mammogram       MRI     Ultrasound      n_t
#001-101 0       65      Yes     36      49      34      7       1       B       3 capsules of 200 mg Ribociclib dai
#001-101 31      55      Yes     NA      NA      NA      7       1       B       3 capsules of 200 mg Ribociclib dai
def read_clinical_all(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        #[ARM, Response]
        data[file_df['Patient.Study.ID'][i]] = [file_df['Response'][i]]
    return data


#Sample  Patient.Study.ID        Arm     ARM     Treatment       Ribo    TreatLab        Class   Response
#FEL001P120      001-120 B       B       letrozole + ribo        1       Letrozole + Ribociclib  PD      Non-responder
#FEL002P106      001-106 C       C       letrozole + ribo        1       Letrozole + Ribociclib  PR      Responder
#FEL003P112      001-112 B       B       letrozole + ribo        1       Letrozole + Ribociclib  SD      Non-responder
def read_clinical(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        #[ARM, Response]
        data[file_df['Patient.Study.ID'][i]] = [file_df['Sample'][i]]
    return data

#Patient.Study.ID        dynamic_class   dynamic_class2
#001-101 Stable disease  Stable disease
#001-101 Stable disease  Stable disease
#001-101 Stable disease  Stable disease
def read_clinical_5groups(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        #[ARM, Response]
        data[file_df['Patient.Study.ID'][i]] = [file_df['dynamic_class'][i], file_df['dynamic_class2'][i]]
    return data

#Sample  Patient.Study.ID        Arm     ARM     Treatment       Ribo    TreatLab        Class   Response
#FEL001P120      001-120 B       B       letrozole + ribo        1       Letrozole + Ribociclib  PD      Non-responder
#FEL002P106      001-106 C       C       letrozole + ribo        1       Letrozole + Ribociclib  PR      Responder
#FEL003P112      001-112 B       B       letrozole + ribo        1       Letrozole + Ribociclib  SD      Non-responder
def add_ki67_to_clinical(infile, ki67_group, v3_group):
    data = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=0)
    Ki67_response = []
    v3_response   = []
    v3_response2  = []
    orig_ident    = []
    for i in range(file_df.shape[0]):
        #[ARM, Response]
        if ki67_group.has_key(file_df['Patient.Study.ID'][i]):
            Ki67_response.append(ki67_group[file_df['Patient.Study.ID'][i]])
        else:
            Ki67_response.append("NoData")
        if v3_group.has_key(file_df['Patient.Study.ID'][i]):
            v3_response.append(v3_group[file_df['Patient.Study.ID'][i]][0])
            v3_response2.append(v3_group[file_df['Patient.Study.ID'][i]][1])
        else:
            v3_response.append("NoData")
            v3_response2.append("NoData")
        orig_ident.append(file_df['Sample'][i][:6])

    file_df['Ki67_Response'] = Ki67_response
    file_df['Response_v2'] = file_df['Response']
    file_df['Response_v3'] = v3_response 
    file_df['Response_v3b'] = v3_response2
    file_df['orig.ident'] = orig_ident
    file_df.drop('Response', axis=1, inplace=True)
    header = ["Sample", "orig.ident", "Patient.Study.ID", "ARM", "Treatment", "Ribo", "Ki67_Response", "Response_v2", "Response_v3", "Response_v3b"]
    file_df= file_df[header]
    outfile = '%s.Ki67_Response.txt' %(infile)
    file_df.to_csv(outfile, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ki67')
    parser.add_argument('--clinical')
    parser.add_argument('--clinical_all')
    parser.add_argument('--clinical_5groups')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.ki67) > 0 and len(args.clinical) > 0
    except:
        usage()
        sys.exit(2)

    # Jason's response v2
    patient_response = read_clinical_all(args.clinical_all)
 
    # jason's response v3 5 groups
    patient_response_v3 = read_clinical_5groups(args.clinical_5groups) 

    # Only FELID -> trial ID
    patient_id = read_clinical(args.clinical)

    # classify ki67 groups and merged with jason's v2
    ki67_group, response_v3 = classify_ki67_group(args.ki67, patient_response, patient_response_v3)
   
    # add ki67 groups to clinical meta
    add_ki67_to_clinical(args.clinical, ki67_group, response_v3)
    


if __name__ == '__main__':
    main()

