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
python scRNA_Seurat_05individual_from_individual_merge_matrix.py --list FEL011027_responder_vs_non_responder.list

--list: list of patients that need to be merged

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

def read_file_pd(infile):
    file_df = pd.read_csv(infile, sep="\t", header=0)
    return file_df

def get_file_df(patients, filetype):
    df_list   = []
    for p in patients:
        filename = '%s%s' %(p, filetype)
        if os.path.exists(filename):
            df_file = read_file_pd(filename)
            df_list.append(df_file)
    return df_list

#FEL027_M_TATATCCTCAAGGTGG
def read_patient_list(infile):
    file_df  = pd.read_table(infile, header=None, names=['Patient'])
    patients = list(file_df['Patient'])
    print patients
    return patients


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--list')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.list) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = "test_merge"

    #file_list = glob.glob('%s/*.ssGSEA.scores.txt' %(args.input))
    #df_list   = []
    #for f in sorted(file_list):
    #    df_file = read_file_pd(f)
    #    df_list.append(df_file)

    patients = read_patient_list(args.list)

    #filetypes = ['_gene_symbols.genefu_filter.sample_sum.txt', '_gene_symbols.genefu.sample_sum.txt', '_gene_symbols.genefu_filter.sample_pc.txt', '_gene_symbols.genefu.sample_pc.txt']
    filetypes = ['_gene_symbols.genefu.sample_sum.txt', '_gene_symbols.genefu.sample_pc.txt']
    for ft in filetypes:
        outfile = '%s.%s' %(args.list, ft)
        if re.compile(r'.list$').search(args.list):
            outfile = re.sub(r'.list', r'%s' %(ft), args.list)
        df_list = get_file_df(patients, ft)
        df_merged = pd.concat(df_list)
        df_merged = df_merged[["Sample", "ER-/HER2-", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif", "HER2+", "Total"]]
        df_merged.to_csv(outfile, sep="\t", index=False)
        outxlsx = '%s.xlsx' %(outfile)
        s = re.compile(r'.txt$')
        if s.search(outfile):
            outxlsx = re.sub(r'.txt', r'.xlsx', outfile)
        df_merged.to_excel(outxlsx, index=False)


    #filetypes = ['_10x_gene_symbols.scaled.counts.txt']
    #for ft in filetypes:
    #    outfile = '%s.%s' %(args.list, ft)
    #    if re.compile(r'.list$').search(args.list):
    #        outfile = re.sub(r'.list', r'%s' %(ft), args.list)
    #    df_list = get_file_df(patients, ft)
    #    df_merged = mergefiles(df_list, on="Gene.ID") 
    #    df_merged.to_csv(outfile, sep="\t", index=False)

    #merge
    #outfile   = '%s.ssGSEA.scores.txt' %(args.output)
    #df_merged = mergefiles(df_list, on="Gene Set")
    #df_merged.to_csv(outfile, sep="\t", index=False)

if __name__ == '__main__':
    main()

