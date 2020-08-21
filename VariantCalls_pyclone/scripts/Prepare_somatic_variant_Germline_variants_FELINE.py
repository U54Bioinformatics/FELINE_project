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
def extract_germline_variants(infile, outfile_variant, cutoff, germline_col):
    germline_col = int(germline_col)
    data = defaultdict(str)
    linen = 0
    header= ''
    chro = map(str, range(1,23))
    chro.append("X")
    chro.append("Y")
    ofile = open(outfile_variant, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            linen += 1
            if linen >= 4:
                #print line
                unit = re.split(r'\t', line)
                #print unit[15:18]
                germline_gt = re.split(r'/', unit[germline_col])
                #60/3/0.046, germline has more than 1 reads
                #60/3/0.046, germline VAF > 1%
                if not len(germline_gt) == 3:
                    # some samples do not have sequence coverage at those sites
                    continue
                #print "pass"
                if float(germline_gt[2]) >= float(cutoff):
                #if float(germline_gt[1]) >= int(cutoff):
                    print >> ofile, '%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], str(int(unit[1])+1), header, unit[germline_col])
                elif not unit[0] in chro:
                    print >> ofile, '%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], str(int(unit[1])+1), 'Nonchr', 'Nonchr')
            elif linen == 3:
                unit = re.split(r'\t', line)
                header = unit[germline_col]
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
    parser.add_argument('--cut')
    parser.add_argument('--germline_col')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0 and len(args.germline_col) > 0
    except:
        usage()
        sys.exit(2)

    # cutoff for minimum read
    # cutoff for max VAF
    if not args.cut: 
        args.cut = 0.01

    # output filename: orig header, raw/unfilter, filter, vcf
    outfile_variant    = '%s_Germline_variant.bed' %(args.input)
    #outfile_raw     = '%s.table.raw.txt' %(args.input)
    #outfile_filter  = '%s.table.filter.txt' %(args.input)
    #outfile_vcf     = '%s.vcf' %(args.input)
    s = re.compile(r'.txt$')
    if s.search(args.input):
        outfile_variant= re.sub(r'.txt', r'.table.germline.txt', args.input)
    #    outfile_raw = re.sub(r'.txt', r'.table.raw.txt', args.input)
    #    outfile_filter = re.sub(r'.txt', r'.table.filter.txt', args.input)
    #    outfile_vcf = re.sub(r'.txt', r'.vcf', args.input)

    # test if VCF header exists and make a copy as header of output VCF file
    #if not os.path.exists(args.vcf_header):
    #    print 'VCF header file does not exist: %s' %(args.vcf_header)
    #else:
    #    os.system('cp %s %s' %(args.vcf_header, outfile_vcf))
    
    # filter variants and output to a VCF file
    extract_germline_variants(args.input, outfile_variant, args.cut, args.germline_col)

if __name__ == '__main__':
    main()

