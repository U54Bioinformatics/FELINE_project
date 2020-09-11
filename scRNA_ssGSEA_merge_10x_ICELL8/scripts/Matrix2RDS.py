#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import re
import os
import argparse
import pyreadr

def usage():
    test="name"
    message='''
python Infercnv_pipeline_individual_from_individual_split_infercnv_input.py --anno FEL025P142_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt

Split large sample into ~5000 cells chunks. Run infercnv on these chunks.

    '''
    print(message)


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

def convert_count_matrix2RDS(infile, outfile):
    file_df = pd.read_csv(infile, sep="\t", header=0)
    pyreadr.write_rds(outfile, file_df)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input')
    parser.add_argument('--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        if args.input.endswith('.txt'): 
            args.output = re.sub(r'.txt', '', args.input)
            args.output = '%s.RDS' %(args.output)

    print(args.input)
    print(args.output)
    convert_count_matrix2RDS(args.input, args.output)

if __name__ == '__main__':
    main()

