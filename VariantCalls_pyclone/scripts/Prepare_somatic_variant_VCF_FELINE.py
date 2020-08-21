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
Extract somatic mutations from BETSY variant calls and save filtered variants in VCF

python Variant_calls_filtered_vcf.py --input 97OVCZ_variant_calls

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

def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data

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
#Output: bed
#1       2583843 2583844 TTC34_1_2583843 FEL013_G
#1       2585401 2585402 TTC34_1_2585401 FEL013_G
#1       2585503 2585504 TTC34_1_2585503 FEL013_G
#1       6027167 6027168 NPHP4_1_6027167 FEL013_G
def variant_vcf(infile, outfile_vcf, outfile_tab):
    data = defaultdict(str)
    linen = 0
    header= ''
    chro = map(str, range(1,23))
    chro.append("X")
    chro.append("Y")
    ofile = open(outfile_vcf, 'a')
    ofile_tab = open(outfile_tab, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            linen += 1
            if linen >= 4:
                #print line
                unit = re.split(r'\t', line)
                #print unit[15:18]
                sample1_gt = convert_vcf_gt(unit[15])
                sample2_gt = convert_vcf_gt(unit[16]) 
                #60/3/0.046, germline has more than 1 reads
                chrom = 'chr%s' %(unit[0])
                pos   = unit[1]
                ref   = unit[2]
                alt   = unit[3]
                if len(ref) == 1 and len(alt) == 1:
                    print >> ofile, '%s\t%s\t.\t%s\t%s\t200\tPASS\tAB=0\tGT:DP:AD:AF\t%s\t%s' %(chrom, pos, ref, alt, sample1_gt, sample2_gt)
                    print >> ofile_tab, '%s\t%s\t.\t%s\t%s\t200\tPASS\tAB=0\tGT:DP:AD:AF\t%s\t%s' %(chrom, pos, ref, alt, sample1_gt, sample2_gt)
            elif linen == 3:
                unit = re.split(r'\t', line)
                header = [unit[15], unit[16]]
                print >> ofile, '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\t%s' %(unit[15], unit[16])
                print >> ofile_tab, 'CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\t%s' %(unit[15], unit[16])
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--vcf_header')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.vcf_header:
        args.vcf_header = 'vcf_header.vcf'


    # output filename: orig header, raw/unfilter, filter, vcf
    outfile_vcf     = '%s.vcf' %(args.input)
    outfile_tab     = '%s.tab' %(args.input)
    s = re.compile(r'.txt$')
    if s.search(args.input):
        outfile_vcf = re.sub(r'.txt', r'.vcf', args.input)
        outfile_tab = re.sub(r'.txt', r'.tab', args.input)

    # test if VCF header exists and make a copy as header of output VCF file
    if not os.path.exists(args.vcf_header):
        print 'VCF header file does not exist: %s' %(args.vcf_header)
    else:
        os.system('cp %s %s' %(args.vcf_header, outfile_vcf))
    
    # output a VCF file
    variant_vcf(args.input, outfile_vcf, outfile_tab)

if __name__ == '__main__':
    main()

