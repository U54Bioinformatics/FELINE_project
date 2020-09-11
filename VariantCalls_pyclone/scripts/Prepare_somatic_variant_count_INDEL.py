#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import pandas as pd

def usage():
    test="name"
    message='''
Count variant numbers of somatic mutations from BETSY variant calls

python Prepare_somatic_variant_count_variant.py --input FELINE_patients.list 

    '''
    print message


def mergefiles(dfs=[], on=''):
    """Merge a list of files based on one column"""
    if len(dfs) == 1:
         return "List only have one element."

    elif len(dfs) == 2:
        df1 = dfs[0]
        df2 = dfs[1]
        df = df1.merge(df2, on=on)
        return df

    # Merge the first and second datafranes into new dataframe
    df1 = dfs[0]
    df2 = dfs[1]
    df = dfs[0].merge(dfs[1], on=on)

    # Create new list with merged dataframe
    dfl = []
    dfl.append(df)

    # Join lists
    dfl = dfl + dfs[2:] 
    dfm = mergefiles(dfl, on)
    return dfm


#FEL013
#FEL015
#FEL019
def read_patient_list(infile):
    patients = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                patients.append(unit[0])
    return patients

#convert besty gt to vcf gt
#60/8/0.118
#GT:DP:AD:AF,
def convert_vcf_gt(betsy_gt):
    gt = re.split(r'/', betsy_gt)
    sample_GT = ''
    if betsy_gt == '':
        return '0/0:0:0:0'
    if int(gt[0]) > 0 and int(gt[1]) > 0:
        sample_GT = '0/1'
    elif int(gt[0]) > 0:
        sample_GT = '0/0'
    else:
        sample_GT = '1/1'    

    sample_DP = str(int(gt[0]) + int(gt[1]))
    sample_AD = '%s,%s' %(gt[0], gt[1])
    sample_AF = gt[2]
    return '%s:%s:%s:%s' %(sample_GT, sample_DP, sample_AD, sample_AF)

#Extract mutations that are present in germline or mutations that need to be removed (not on chr1 to 22, X, Y)
#Input: BESTY flat table
#                                Cancer Genes                                                                                    Coverage                        Annovar
#
#Chrom   Pos     Ref     Alt     COSMIC  Stratton        TARGET  Vogelstein      Campbell        Thunderbolts    Pancreatic Cancer (Scarlett)    Achilles (Cheung)       Pancreatic
#1       809319  C       T                                                                                               60/8/0.118      29/0/0.000      57/7/0.109      ncRNA_intronic  F
#1       1904429 G       A                                                                                               34/0/0.000      58/0/0.000      42/8/0.160      exonic  CFAP74
#1       2583843 A       T                                                                                               39/2/0.049      60/3/0.046      25/7/0.219      intronic
def INDEL_count(infile):
    count_total = 0
    linen = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            linen += 1
            if linen >= 4:
                unit = re.split(r'\t', line)
                try:
                    if len(unit[2]) > 1 or len(unit[3]) > 1:
                        count_total += 1
                except:
                    continue
    return count_total

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input')
    parser.add_argument('--variant_folder')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0 & len(args.variant_folder) > 0
    except:
        usage()
        sys.exit(2)

    patients = read_patient_list(args.input)
    folder   = args.variant_folder 

    print 'Patient\tINDEL'
    for p in patients:
        svm_file = '%s/%s_variant_calls_04_purity40_totalRD20_minread5_minVAF0.05.txt' %(folder, p)
        indel    = INDEL_count(svm_file)    
        print '%s\t%s' %(p, indel)
if __name__ == '__main__':
    main()

