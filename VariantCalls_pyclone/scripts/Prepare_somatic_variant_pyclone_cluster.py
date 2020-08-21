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

#sample_id       cluster_id      size    mean    std
#FEL045_S        0       1189    0.499974927780264       0.0005854684879395742
#FEL045_M        0       1189    0.4902703097562532      0.0016217930188956308
#FEL045_S        1       70      0.46851300759812015     0.012818326299048012
#FEL045_M        1       70      2.2442932943855042e-13  4.7373972752294485e-08
def create_pyclone_final_cluster(infile):
    samples        = defaultdict(lambda : int())
    cluster_info   = defaultdict(lambda : defaultdict(lambda: str()))
    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        #print file_df[0][i]
        file_df['size'][i]
        cluster_info[file_df['cluster_id'][i]]["size"]                   = file_df['size'][i]
        cluster_info[file_df['cluster_id'][i]][file_df['sample_id'][i]]  = file_df['mean'][i]
        samples[file_df['sample_id'][i]] = 1 
    ## write final cluster file
    outfile = re.sub(r'.txt', r'.final.txt', infile)
    ofile   = open(outfile, 'w')
    sample_names = sorted(samples.keys())
    print >> ofile, "ClusterID\tCluster\tNum Mutations\t%s" %('\t'.join(sample_names))
    for cluster in sorted(cluster_info.keys()):
        if int(cluster_info[cluster]["size"]) > 1:
            new_list = ["1", str(cluster), str(cluster_info[cluster]["size"])]
            for sample in sample_names:
                new_list.append(str(cluster_info[cluster][sample]))
            print >> ofile, '\t'.join(new_list)
    ofile.close()

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

    create_pyclone_final_cluster(args.input)

if __name__ == '__main__':
    main()

