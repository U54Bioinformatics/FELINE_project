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
    header  = ['Cell.ID', 'HMM_infer']
    file_df = file_df[header]
    return file_df

def get_file_df(patients, filetype):
    df_list   = []
    for p in patients:
        filename = '%s%s' %(p, filetype)
        if os.path.exists(filename):
            df_file = read_file_pd(filename)
            df_list.append(df_file)
    return df_list

#FEL011P138
#FEL012P101
#FEL013P102
def read_patient_list(infile):
    file_df  = pd.read_table(infile, header=None, names=['Patient'])
    patients = file_df['Patient']
    print patients
    return patients

def read_patient_list_timepoint(infile):
    file_df  = pd.read_table(infile, header=None, names=['Patient'])
    patients = []
    for p in file_df['Patient']:
        patients.append('%s_S' %(p))
        patients.append('%s_M' %(p))
        patients.append('%s_E' %(p))
    print patients
    return patients


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--list')
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

    #FEL014P103_10x_cell_metadata.scater_QC.txt
    filetypes = ['_HMM_infer.subclone.anno.txt']
    for ft in filetypes:
        outfile = '%s%s' %(args.output, ft)
        df_list = get_file_df(patients, ft)
        df_merged = pd.concat(df_list)
        df_merged.to_csv(outfile, sep="\t", index=False)


    #filetypes = ['_10x_gene_symbols.scaled.counts.normalized_to_control.txt', '_10x_gene_symbols.scaled.Cancer_cells.mean.txt', '_10x_gene_symbols.scaled.Endothelial_cells.mean.txt', '_10x_gene_symbols.scaled.Fibroblasts.mean.txt', '_10x_gene_symbols.scaled.Macrophages.mean.txt']
    #filetypes = ['_gene_symbols.normalized.counts.txt', '_gene_symbols.scaled.counts.txt']
    #for ft in filetypes:
    #    print ft
    #    outfile = '%s%s' %(args.prefix, ft)
    #    df_list = get_file_df(patients, ft)
    #    df_merged = mergefiles(df_list, on="Gene.ID") 
    #    df_merged.to_csv(outfile, sep="\t", index=False)

    #merge
    #outfile   = '%s.ssGSEA.scores.txt' %(args.output)
    #df_merged = mergefiles(df_list, on="Gene Set")
    #df_merged.to_csv(outfile, sep="\t", index=False)

if __name__ == '__main__':
    main()

