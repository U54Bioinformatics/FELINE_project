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
Extract dbsnp info for variants

python Prepare_somatic_variant_common_variant_FELINE.py --input FEL046_variant_calls_Germline_free --dbsnp_col 26
python Prepare_somatic_variant_common_variant_FELINE.py --input FEL046_variant_calls_04_purity40_totalRD15_minread3_minVAF0.05.txt --dbsnp_col 25

--dbsnp_col: 26 for three samples, 25 for two samples

awk -F"\t" '{print NF}' FEL046_variant_calls_04_purity40_totalRD15_minread3_minVAF0.05.common_snp_info.txt | less -S
awk -F"\t" '$11>0' FEL046_variant_calls_04_purity40_totalRD15_minread3_minVAF0.05.common_snp_info.txt | less -S


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
def extract_dbsnp_info(infile, outfile, dbsnp_col):
    dbsnp_col = int(dbsnp_col)
    freq_col  = dbsnp_col - 1
    linen = 0
    header= ''
    chro = map(str, range(1,23))
    chro.append("X")
    chro.append("Y")
    ofile = open(outfile, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            linen += 1
            if linen >= 4:
                #print line
                unit = re.split(r'\t', line)
                if unit[0] in chro:
                    sample_n     = dbsnp_col - 23
                    out_list = []
                    # variant position
                    out_list.extend(unit[0:4])
                    # variant genotype
                    out_list.extend(unit[15:15+sample_n])
                    # annotation
                    out_list.extend(unit[dbsnp_col-8:dbsnp_col-4])
                    # dbsnp
                    out_list.extend(unit[freq_col:dbsnp_col+1])
                    print >> ofile, '\t'.join(out_list)
            elif linen == 3:
                unit = re.split(r'\t', line)
                sample_n = dbsnp_col - 23
                out_list = []
                out_list.extend(unit[0:4])
                out_list.extend(unit[15:15+sample_n])
                out_list.extend(unit[dbsnp_col-8:dbsnp_col-4])
                out_list.extend(unit[freq_col:dbsnp_col+1])
                print >> ofile, '\t'.join(out_list)
    ofile.close()

#def extract_germline_variants(infile, outfile_variant):
#    file_df = pd.read_csv(infile, sep="\t", header=[0,1,2])
#    print file_df.columns
    #file_df_subset = file_df.iloc[0:file_df.shape[0],[0,1,2,3,4,5,7,42,43,44,45,46,76,77,78,79,80,81]]
    # output a table of raw/unfiltered variants with original header
    #file_df_subset.to_csv(outfile_orig, sep="\t", index=False)
    # this new header need to change for each project. There may have difference in number of samples. So column numbers will be different too.
    #file_df_subset.columns = ["Chrom","Pos","Ref","Alt","Func.refGene","Gene.refGene","ExonicFunc.refGene","Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID","Num_Caller_AB11A","Num_Caller_AB11J","Num_Caller_Germline","Coverage_AB11A","Coverage_AB11J","Coverage_Germline"]
    # output a table of raw/unfiltered variants 
    #file_df_subset.to_csv(outfile_raw, sep="\t", index=False)

    # filter variants: 
    # 1. total read coverage >= 20 in all samples; 
    # 2. alternative alleles read coverate >= 5 in tumors and <= 1 in normal;
    # 3. allele frequency > 0.05?
    # 4. not in chrom 1 to X; not SNPs?;
    #chrs  = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--dbsnp_col')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.dbsnp_col:
        args.dbsnp_col = 23

    # output filename: orig header, raw/unfilter, filter, vcf
    outfile_variant    = '%s.common_snp_info.txt' %(args.input)
    s = re.compile(r'.txt$')
    if s.search(args.input):
        outfile_variant= re.sub(r'.txt', r'.common_snp_info.txt', args.input)

    # extract esp6500siv2_all and snp138 column to check if germline free still have common snp 
    extract_dbsnp_info(args.input, outfile_variant, args.dbsnp_col)

if __name__ == '__main__':
    main()

