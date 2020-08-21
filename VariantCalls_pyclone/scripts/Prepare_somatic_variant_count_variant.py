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
def variant_count(infile, col):
    count_total = 0
    count_exonic=0
    linen = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            linen += 1
            if linen >= 4:
                #print line
                count_total += 1
                unit = re.split(r'\t', line)
                try:
                    #print unit[16], unit[17], unit[18]
                    variant_pos = unit[int(col)]
                    if variant_pos == 'exonic':
                        count_exonic += 1
                except:
                    #print line
                    #sys.exit(2)
                    continue
    return [str(count_total), str(count_exonic)]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--RD')
    parser.add_argument('--minread')
    parser.add_argument('--minVAF')
    parser.add_argument('--mincaller')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    patients = read_patient_list(args.input)

    min_total_depth = 20
    min_alt_depth   = 5
    min_vaf         = 0.05
    min_caller      = 1
    if args.RD:
        min_total_depth = int(args.RD)
    if args.minread:
        min_alt_depth   = int(args.minread)
    if args.minVAF:
        min_vaf         = float(args.minVAF)
    if args.mincaller:
        min_caller      = int(args.mincaller)

    print 'Patient\tTotal\tTotal_exonic\tGermline_Free\tGermline_Free_exonic\tRD\tRD_exonic\tminread\tminread_exonic\tminVAF\tminVAF_exonic\tmincaller\tmincaller_exonic'
    # count variant for each patients
    for p in patients:
        variant_total   = '%s_variant_calls' %(p)
        variant_gl_free = '%s_variant_calls_01_purity40.txt' %(p)
        variant_gl_free_tdepth = '%s_variant_calls_02_purity40_totalRD%s.txt' %(p, min_total_depth)
        variant_gl_free_adepth = '%s_variant_calls_03_purity40_totalRD%s_minread%s.txt' %(p, min_total_depth, min_alt_depth)
        variant_gl_free_adepth_vaf = '%s_variant_calls_04_purity40_totalRD%s_minread%s_minVAF%s.txt' %(p, min_total_depth, min_alt_depth, min_vaf)
        variant_gl_free_adepth_vaf_caller = '%s_variant_calls_05_purity40_totalRD%s_minread%s_minVAF%s_mincaller%s.txt' %(p, min_total_depth, min_alt_depth, min_vaf, min_caller)
        # count raw variants with gerlmine, col=18 is variant position (exonic or not)
        count_sum = []
        count_total, count_exonic = variant_count(variant_total, 18)
        count_sum.append(count_total)
        count_sum.append(count_exonic)
        # count gerlmine free variants, col=17 is variant position (exonic or not)
        for variant in [variant_gl_free, variant_gl_free_tdepth, variant_gl_free_adepth, variant_gl_free_adepth_vaf, variant_gl_free_adepth_vaf_caller]:
            count_total, count_exonic = variant_count(variant, 17)
            count_sum.append(count_total)
            count_sum.append(count_exonic)
        print '%s\t%s' %(p, '\t'.join(count_sum))
        

if __name__ == '__main__':
    main()

