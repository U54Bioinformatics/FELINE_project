#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import re
import os
import argparse
import glob

def usage():
    test="name"
    message='''
python Prepare_subtype_per_sample.py --subtype LumA --patient_data FEL011040_patients_gene_symbols.genefu_filter.sample_sum.txt --patient_responder ../Celltype/FEL011_039.responder.txt

Input:
--patient_data:
Sample  Basal   Her2    LumA    LumB    Normal  Total
FEL011P138_M            1.0     82      16.0    1.0     100
FEL012P101_E    8.0     22.0    387     18.0    27.0    462
FEL012P101_M    1.0     5.0     111     14.0    1.0     132
FEL013P102_E    21.0    129.0   439     44.0    84.0    717
FEL013P102_M    11.0    126.0   448     27.0    40.0    652
FEL013P102_S    132.0   16.0    491     23.0    34.0    696
FEL014P103_E    0.0     1.0     4       1.0     0.0     6
FEL014P103_S    5.0     18.0    104     66.0    5.0     198

--patient_responder:
Patient.Study.ID        Sample  MRI_response    Clinical_response       prop_change     Arm     ARM     Ribo    Treatment       TreatLab        Class   Response        Burden_t0
001-101 FEL012P101      SD      NC      0.865028111458766       B       B       1       letrozole + ribo        Letrozole + Ribociclib  SD      Non-responder   49.0050693012157
001-102 FEL013P102      SD      CR      0.668652354145204       C       C       1       letrozole + ribo        Letrozole + Ribociclib  PR      Responder       18.0537794823837
001-103 FEL014P103      SD      PR      0.887085404117355       A1      A       0       letrozole       Letrozole       SD      Non-responder   68.6648264425512        ControlArm
001-104 FEL010P104      SD      CR      0.722067100177336       B       B       1       letrozole + ribo        Letrozole + Ribociclib  SD      Non-responder   24.9986814409199
001-105 FEL015P105      SD      CR      0.843149611510377       A2      A       0       letrozole       Letrozole       SD      Non-responder   13.0159549065415        ControlArm



Output:
Patient Arm     Response        Macrophages_S   Macrophages_M   Macrophages_E   Change_M_S      Change_E_S
FEL011P138      B       nan     2.67295597484   2.8818443804    0       0.208888405561  NA
FEL012P101      B       NC      3.25264750378   4.03022670025   0.39481105471   0.77757919647   -2.85783644907
FEL013P102      C       CR      1.03469263542   2.03449800973   8.94368789106   0.999805374307  7.90899525563
FEL014P103      A       PR      9.3392395718    10.6522759268   6.08465608466   1.313036355     -3.25458348714


    '''
    print message

#Sample  Basal   Her2    LumA    LumB    Normal  Total
#FEL011P138_M            1.0     82.0    16.0    1.0     100
#FEL012P101_E    1.7316017316    4.7619047619    83.7662337662   3.8961038961    5.84415584416   462
#FEL012P101_M    0.757575757576  3.78787878788   84.0909090909   10.6060606061   0.757575757576  132
def read_patient_subtype(infile):
    data = defaultdict(lambda : defaultdict(lambda : list()))
    file_df = pd.read_table(infile, header=0)
    print "Raw data:"
    print file_df
    file_df_drop   = file_df[file_df['Total'] < 10]
    file_df_drop   = file_df_drop.reset_index(drop=True)
    file_df_filter = file_df[file_df['Total'] >= 10]
    file_df_filter = file_df_filter.reset_index(drop=True)
    print "After filter by 50 cells:"
    print file_df_filter
    #remove total and calculate percent
    del file_df_filter['Total']
    for i in range(file_df_filter.shape[0]):
        #print file_df_filter.loc[i, 'Sample']
        unit = re.split(r'_', file_df_filter.loc[i, 'Sample'])
        patient   = unit[0]
        timepoint = unit[1]
        print patient, timepoint
        #calculate percent
        percent = map(str, file_df_filter.iloc[i,1:])
        print percent
        #percent = file_df_filter.iloc[i,2:]/file_df_filter.iloc[i,2:].sum() 
        #percent = percent.apply(lambda x: x/x.sum())
        #print percent
        data[patient][timepoint] = percent
    # dropped patient
    for i in range(file_df_drop.shape[0]):
        print i
        unit = re.split(r'_', file_df_drop.loc[i, 'Sample'])
        patient   = unit[0]
        timepoint = unit[1]
        print patient, timepoint
        data[patient][timepoint] = ['nan', 'nan', 'nan', 'nan']
    return data

def patient_subtype_matrix(subtype_data, output): 
    ofile = open(output, 'w') 
    timepoint = ['S', 'M', 'E']
    subtype   = ['Basal', 'LumB', 'LumA', 'Her2+']
    print >> ofile, 'Sample\t%s\t%s\t%s' %('\t'.join(['S_%s' %(x) for x in subtype]), '\t'.join(['M_%s' %(x) for x in subtype]), '\t'.join(['E_%s' %(x) for x in subtype]))
    for p in sorted(subtype_data.keys()):
        values = [p]
        for t in timepoint:
            if subtype_data[p].has_key(t):
                values.extend(subtype_data[p][t])
            else:
                values.extend(['nan', 'nan', 'nan', 'nan'])
        print >> ofile, '%s' %('\t'.join(values))
    ofile.close()

#Patient.Study.ID        Sample  MRI_response
def read_patient_responder(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        patient  = file_df['Sample'][i]
        responder= file_df['Response'][i]
        treat    = file_df['TreatCode'][i]
        data[patient] = [responder, treat]
    return data
  
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subtype')
    parser.add_argument('--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.subtype) > 0
    except:
        usage()
        sys.exit(2)
   
    if not args.output:
        args.output = 'FEL011043.subtype_heatmap_matrix.txt'
 
    subtype_data = read_patient_subtype(args.subtype)
    patient_subtype_matrix(subtype_data, args.output)

if __name__ == '__main__':
    main()

