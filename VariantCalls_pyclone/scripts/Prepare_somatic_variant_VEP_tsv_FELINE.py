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

#'Chromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tt_depth\tt_ref_count\tt_alt_count\tn_depth\tn_ref_count\tn_alt_count'
def convert_tsv_line(chro, pos, ref, alt, gt_G, gt_S, t_barcode, n_barcode):
    if len(gt_S) == 1:
        gt_S = [0, 0]
    if len(gt_G) == 1:
        gt_G = [0, 0]

    t_allele1 = ref if int(gt_S[0]) > 0 else '-'
    t_allele2 = alt if int(gt_S[1]) > 0 else '-'
    tt_depth  = int(gt_S[0]) + int(gt_S[1])
    tt_ref    = int(gt_S[0])
    tt_alt    = int(gt_S[1])
    tn_depth  = int(gt_G[0]) + int(gt_G[1])
    tn_ref    = int(gt_G[0])
    tn_alt    = int(gt_G[1])
    n_allele1 = ref if int(gt_G[0]) > 0 else '-'
    n_allele2 = alt if int(gt_G[1]) > 0 else '-'
    return [chro, pos, ref, t_allele1, t_allele2, t_barcode, n_barcode, n_allele1, n_allele2, tt_depth, tt_ref, tt_alt, tn_depth, tn_ref, tn_alt] 

#Extract mutations that are present in germline or mutations that need to be removed (not on chr1 to 22, X, Y)
#Input: BESTY flat table
#                                Cancer Genes                                                                                    Coverage                        Annovar
#
#Chrom   Pos     Ref     Alt     COSMIC  Stratton        TARGET  Vogelstein      Campbell        Thunderbolts    Pancreatic Cancer (Scarlett)    Achilles (Cheung)       Pancreatic 
#1       809319  C       T                                                                                               60/8/0.118      29/0/0.000      57/7/0.109      ncRNA_intronic  F
#1       1904429 G       A                                                                                               34/0/0.000      58/0/0.000      42/8/0.160      exonic  CFAP74
#1       2583843 A       T                                                                                               39/2/0.049      60/3/0.046      25/7/0.219      intronic        
#Output: vep.tsv
#Chromosome	Start_Position	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	Tumor_Sample#_Barcode	Matched_Norm_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2	t_depth	t_ref_countt_alt_count	n_depth	n_ref_count	n_alt_count
#1	11290179	A	A	-	TUMOR	NORMAL	A	A	20	10	10	35	30	5
def extract_variants_to_vep_tsv(snp_raw, snp_filtered, outfile_variant, patient, germline_col):
    germline_col = int(germline_col)
    data = defaultdict(str)
    linen = 0
    header= ''
    ofile = open(outfile_variant, 'w')
    print >> ofile, 'Chromosome\tStart_Position\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tt_depth\tt_ref_count\tt_alt_count\tn_depth\tn_ref_count\tn_alt_count'
    with open (snp_raw, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            linen += 1
            if linen >= 4:
                unit = re.split(r'\t', line)
                snp  = '%s_%s'%(unit[0], unit[1])
                ref  = unit[2]
                alt  = unit[3]
                if snp_filtered.has_key(snp):
                    gt_E = []
                    gt_G = []
                    gt_S = []
                    sample_E = '%s_E' %(patient)
                    sample_G = '%s_G' %(patient)
                    sample_S = '%s_S' %(patient)
                    if germline_col == 16:
                        sample_E = '%s_E' %(patient)
                        gt_E = re.split(r'/', unit[15])
                        gt_G = re.split(r'/', unit[16])
                        gt_S = re.split(r'/', unit[17])
                    elif germline_col == 15:
                        sample_E = '%s_M' %(patient)
                        gt_E = re.split(r'/', unit[16])
                        gt_G = re.split(r'/', unit[15])
                        gt_S = re.split(r'/', unit[17])
                    tsv_line_S = convert_tsv_line(unit[0], unit[1], ref, alt, gt_G, gt_S, sample_S, sample_G)
                    tsv_line_E = convert_tsv_line(unit[0], unit[1], ref, alt, gt_G, gt_E, sample_E, sample_G)
                    print >> ofile, '\t'.join(map(str, tsv_line_S))
                    print >> ofile, '\t'.join(map(str, tsv_line_E))
    ofile.close()

def extract_filtered_variants(infile):
    data = defaultdict(str)
    linen = 0
    header= ''
    chro = map(str, range(1,23))
    chro.append("X")
    chro.append("Y")
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            linen += 1
            if linen >= 4:
                #print line
                unit = re.split(r'\t', line)
                snp  = '%s_%s' %(unit[0], unit[1])
                data[snp] = 1
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input')
    parser.add_argument('--filtered')
    parser.add_argument('--patient')
    parser.add_argument('--germline_col')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.patient) > 0
    except:
        usage()
        sys.exit(2)

    if not args.germline_col:
        args.germline_col = 16

    if not args.input:
        args.input = "%s_variant_calls" %(args.patient)

    if not args.filtered:
        #args.filtered = "%s_variant_calls_04_purity40_totalRD15_minread3_minVAF0.05.txt" %(args.patient)
        args.filtered = "%s_variant_calls_04_purity40_totalRD20_minread5_minVAF0.05.txt" %(args.patient)

    outfile_variant    = '%s.vep.tsv' %(args.filtered)
    s = re.compile(r'.txt$')
    if s.search(args.filtered):
        outfile_variant= re.sub(r'.txt', r'.vep.tsv', args.filtered)
  
    # read filtered snp 
    snp_filtered = extract_filtered_variants(args.filtered)
 
    # extract info from raw snp file
    extract_variants_to_vep_tsv(args.input, snp_filtered, outfile_variant, args.patient, args.germline_col)

if __name__ == '__main__':
    main()

