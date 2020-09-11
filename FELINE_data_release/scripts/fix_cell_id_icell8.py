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

#Cell.ID	orig.ident	nCount_RNA	nFeature_RNA	Sample	Timepoint	Day
#FEL001_C25287_104812_C00_R10	FEL001	39374.0	2640	FEL001P120	S	0
#FEL001_C25287_104812_C00_R17	FEL001	101195.01	2178	FEL001P120	S	0
#FEL001_C25287_104812_C00_R46	FEL001	600975.99	4808	FEL001P120	S	0
def new_cell_id_icell8(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    outfile = "FEL001010_Cancer_cells_icell8.cell_id_fixed.txt"
    ofile   = open(outfile, 'w')
    for i in range(file_df.shape[0]):
        if file_df['Platform'][i] == "ICELL8":
            unit = re.split(r'_', file_df['Cell.ID'][i])
            cell_id_list = [file_df['orig.ident'][i], file_df['Timepoint'][i]] + unit[1:]
            cell_id_new  = '_'.join(cell_id_list)
            data[file_df['Cell.ID'][i]] = cell_id_new
            print >> ofile, "%s\t%s" %(file_df['Cell.ID'][i], cell_id_new) 
    ofile.close()
    return data

def replace_cell_id_icell8(infile, cell_id):
    data = defaultdict(lambda : str())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    outfile = re.sub(r'.txt', '.id_fixed.txt', infile)
    outfile1= "FEL001010_Cancer_cells_icell8.cell_id_replaced.txt"
    ofile1  = open(outfile1, 'w')
    # replace header
    header_old = file_df.columns
    header_new = []
    for h in header_old:
        h_new = ''
        if h == "Gene.ID":
            header_new.append(h)
            h_new = "Gene.ID"
        else:
            if cell_id.has_key(h):
                header_new.append(cell_id[h])
                h_new = cell_id[h]
            else:
                header_new.append("NA")
                h_new = "NA"
        print >> ofile1, "%s\t%s" %(h, h_new)
    ofile1.close()
    file_df.columns = header_new
    file_df.to_csv(outfile, sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--meta')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0 and len(args.meta) > 0
    except:
        usage()
        sys.exit(2)

    # generated new cell id
    cell_id = new_cell_id_icell8(args.meta)
    # replace cell id from count file
    replace_cell_id_icell8(args.input, cell_id)

if __name__ == '__main__':
    main()

