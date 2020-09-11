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

#Cell.ID orig.ident      nCount_RNA      nFeature_RNA    Sample  Timepoint       Day     Celltype
def reformat_meta(infile):
    file_df  = pd.read_table(infile, header=0)
    new_anno = []
    for i in range(file_df.shape[0]):
        #print file_df[0][i]
        if file_df["Celltype"][i] in ["Adipocytes", "Endothelial cells", "Fibroblasts", "Pericytes"]:
            new_anno.append("Stromal cells")
        elif file_df["Celltype"][i] in ["B cells", "T cells", "Macrophages", "HSC"]:
            new_anno.append("Immune cells")
        else:
            new_anno.append(file_df["Celltype"][i])
    file_df["Celltype_orig"] = file_df["Celltype"]
    file_df["Celltype"] = new_anno
    return file_df

def read_large_matrix(infile):
    # determine and optimize dtype
    # Sample 100 rows of data to determine dtypes.
    file_test = pd.read_csv(infile, sep="\t", header=0, nrows=100)
    float_cols = [c for c in file_test if file_test[c].dtype == "float64"]
    int_cols = [c for c in file_test if file_test[c].dtype == "int64"]
    if float_cols > 0:
        dtype_cols = {c: np.float32 for c in float_cols}
    elif int_cols > 0:
        dtype_cols = {c: np.int32 for c in int_cols}
    file_df = pd.read_csv(infile, sep="\t", header=0, dtype=dtype_cols)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    file_df = reformat_meta(args.input)
    outfile = "FEL011046.cell_metadata.txt"
    file_df.to_csv(outfile, sep="\t", index=False) 

if __name__ == '__main__':
    main()

