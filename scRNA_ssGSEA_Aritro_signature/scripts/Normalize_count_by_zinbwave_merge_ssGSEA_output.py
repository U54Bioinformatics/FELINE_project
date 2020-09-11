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
python Infercnv_pipeline_individual_from_individual_merge_infercnv_obs.py --subsample FEL028P131a,FEL028P131b

Merge infercnv results from FEL028P131a, FEL028P131b into FEL028P131.
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
    #file_df["Gene.ID"] = list(file_df.index.values)
    return file_df

def read_file_list(infile):
    file_df = pd.read_csv(infile, sep="\t", header=None)
    file_df.columns = ['File']
    #file_df["Gene.ID"] = list(file_df.index.values)
    return file_df['File']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subfiles')
    parser.add_argument('--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.subfiles) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = 'zinbwave_merge'

    #parse and prepare input and output filenames
    file_chunks = read_file_list(args.subfiles)

    ofile_log = open('%s.log' %(args.output), 'w')
    #merge observation files
    zbfiles = ['zinbwave.normalized.ssGSEA_scores.txt']
    #obsfiles = ['test.data.txt']
    for zb in zbfiles:
        print >> ofile_log, 'processing %s type files' %(zb)
        zb_df_list = []
        for file_name in file_chunks:
            zb_filename = re.sub(r'txt', r'%s' %(zb), file_name)
            print >> ofile_log, 'reading %s file' %(zb_filename)
            file_df   = read_file_pd(zb_filename)
            print >> ofile_log, 'File dimension %s %s' %(str(file_df.shape[0]), str(file_df.shape[1]))
            zb_df_list.append(file_df)
        zb_merged_file = '%s.%s' %(args.output, zb)
        print >> ofile_log, 'merging into %s' %(zb_merged_file)
        zb_df_merged = mergefiles(zb_df_list, on="Gene Set")
        zb_df_merged.set_index("Gene Set", inplace=True)
        #zb_df_merged.drop("Gene Set", inplace=True)
        print >> ofile_log, 'File dimension %s %s' %(str(zb_df_merged.shape[0]), str(zb_df_merged.shape[1]))
        #write output into file with orignial name
        zb_df_merged.to_csv(zb_merged_file, sep="\t", index=True)


if __name__ == '__main__':
    main()
