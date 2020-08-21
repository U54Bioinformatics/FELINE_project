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
python CircosConf.py --input circos.config --output pipe.conf
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
def read_patient(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = 1
    return data

def read_patient(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = 1
    return data

#Normal  Cancer
#FEL013_G        FEL013_E
#FEL013_G        FEL013_S
#FEL015_G        FEL015_E
#Or
#Filename        Sample  Group   Pair
#SH5647_SA43217_S1_L004_R1_001.fastq     FEL013_S        single_group    1
#SH5647_SA43217_S1_L004_R2_001.fastq     FEL013_S        single_group    2
def subset_table(infile, patient, col, coln, name, tab):
    samples = ['S', 'E', 'G']
    samples = ['%s_%s' %(patient, x) for x in samples]
    samples.append(col)
    #print samples
    ofile = open('%s%s' %(patient, name), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\s+',line)
                if unit[int(coln)] in samples:
                    if int(tab) == 1:
                        print >> ofile, '\t'.join(unit)
                    else:
                        print >> ofile, '     '.join(unit)
    ofile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--patient')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.patient) > 0
    except:
        usage()
        sys.exit(2)

    home='/home/jichen/Projects/Breast/scRNA/FELINE/FELINE_merged_011_040_premrna/VariantCalls_SNPs'
    patients = read_patient(args.patient)
    for p in sorted(patients):
        subset_table('FELINE_sample.txt', p, 'Sample', 1, '_sample.txt', 1)
        subset_table('FELINE_normal_cancer.txt', p, 'Cancer', 1, '_normal_cancer.txt', 1)
        bam_dir = '%s_bam' %(p)
        if not os.path.exists(bam_dir):
            os.mkdir(bam_dir)
        ln_bam = 'ln -s %s/FELINE_bam/%s*.bam %s/' %(home, p, bam_dir) 
        ln_bai = 'ln -s %s/FELINE_bam/%s*.bam.bai %s/' %(home, p, bam_dir)
        os.system(ln_bam)
        os.system(ln_bai)

if __name__ == '__main__':
    main()

