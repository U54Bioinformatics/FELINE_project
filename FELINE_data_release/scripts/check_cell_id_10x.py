#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import pandas
import re
import os
import argparse
import glob

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

#FEL024_S_AACGGGATCATGTCAG  
#FEL024_S_AAAGTGACATGCCGAC   
#FEL024_S_ACTTCGCGTATACCCA
def check_cell_id_10x(infile):
    data = defaultdict(lambda : defaultdict(lambda : int()))
    file_df = pd.read_csv(infile, header=0, sep="\t")
    file_df.columns = ["Cell.ID"]
    for i in range(file_df.shape[0]):
        unit = re.split(r'_', file_df['Cell.ID'][i])
        data[unit[0]][unit[1]] = 1
    for p in sorted(data.keys()):
        for tp in sorted(data[p]):
            print "%s\t%s" %(p, tp)
   
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--platform')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.platform:
        args.platform = '10x'

    if args.platform == '10x':
        check_cell_id_10x(args.input)

if __name__ == '__main__':
    main()

