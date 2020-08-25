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
python Infercnv_pipeline_individual_from_individual_split_infercnv_input.py --anno FEL025P142_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt

Split large sample into ~5000 cells chunks. Run infercnv on these chunks.

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

#FEL025P142_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt
#FEL025P142_10x_gene_symbols.filtered.counts.txt
def split_anno(infile, countfile):
    patient = re.split('_', infile) 
    count   = re.split('_', countfile)
    file_df = pd.read_csv(infile, sep="\t", header=0)
    line_num    = file_df.shape[0]
    line_target = 2000
    line_filen  = (int(line_num)/line_target) + 1
    line_actrual= int(line_num)/line_filen
    suffix       = ['_chunk%s' %(x) for x in range(100)] 
    line_suffix  = suffix[:line_filen]

    print 'Line numbers: %s' %(line_num) 
    print 'Target size: %s' %(line_target)
    print 'Actrual size: %s' %(line_actrual)
    print 'File number: %s' %(line_filen)
    print 'File suffix: %s' %(line_suffix)
    # shuffler row in df
    file_df = file_df.sample(frac=1).reset_index(drop=True)
    # split into chunks
    file_df_split = np.array_split(file_df, line_filen)
    for i in range(len(file_df_split)):
        file_out = '%s%s_%s' %(patient[0], line_suffix[i], '_'.join(patient[1:]))
        count_out= '%s%s_%s' %(patient[0], line_suffix[i], '_'.join(count[1:]))
        print file_out
        print count_out
        file_df_split[i].to_csv(file_out, sep="\t", index=False)        
        #count file
        if os.path.exists(countfile):
            cmd_ln = 'ln -s %s %s' %(countfile, count_out)
            os.system(cmd_ln)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--anno')
    parser.add_argument('--count')
    parser.add_argument('--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.anno) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = "infercnv_split"

    split_anno(args.anno, args.count)

if __name__ == '__main__':
    main()

