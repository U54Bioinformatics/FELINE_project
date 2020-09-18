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
python Compare_cell_type_annotation_ind_vs_merge.py --ind FEL021P118_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt --merge FEL011024_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def clean_name(cell):
    unit = re.split(r'_', cell)
    newname = '_'.join(unit[1:])
    return newname

#Cell.ID orig.ident      nCount_RNA      nFeature_RNA    Sample  Total.Reads     Total.Reads..Thousands. Total.Non.Mitochondria.Reads    Total.Non.Mitochondria.Reads..Thousands.
def read_file_pd(infile, title, clean):
    s = re.compile(r'subclone(\d+)')
    m = s.search(infile)
    subclone = 'HMM_inferk7'
    if m:
        n = m.groups(0)[0]
        subclone = 'HMM_inferk%s' %(n)
    file_df = pd.read_table(infile, header=0)
    #for i in range(file_df.shape[0]):
    #    data[file_df["Cell.ID"][i]] = file_df["Encode_main_cluster"][i]
    #extract on these two columns
    file_df_anno = file_df[["Cell.ID", subclone]]
    #rename column names
    file_df_anno.columns = ["Cell.ID", '%s_%s' %(subclone, title)]
    #remove prefix from names
    if clean:
        file_df_anno["Cell.ID"] = file_df_anno["Cell.ID"].apply(clean_name)
    subclone = '%s_%s' %(subclone, title)
    return file_df_anno, subclone

def read_info(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                unit[0] = re.sub(r'FEL021P118_', r'', unit[0])
                data[unit[0]] = '\t'.join(unit[1:])
    return data


def main():
   
    parser = argparse.ArgumentParser()
    parser.add_argument('--file1')
    parser.add_argument('--file2')
    parser.add_argument('--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.file1) > 0 and len(args.file2) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = "test"

    anno_1_df, title1 = read_file_pd(args.file1, "file1", 0)
    anno_2_df, title2 = read_file_pd(args.file2, "file2", 0)
    anno_df = pd.merge(anno_1_df, anno_2_df, on='Cell.ID', how='inner')
   # print anno_df[1:3,1:3]
    cross = pd.crosstab(anno_df[title1], anno_df[title2])
    cross = cross.T
    #cross.to_csv("%s_compare.csv" %(args.output), sep=",")
    print cross

if __name__ == '__main__':
    main()
