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
python scRNA_add_cell_type_immune_icell8.py --meta FEL001010_icell8.immune.metadata.txt

Annoate immune cells in FELINE icell8 data using SingleR annoateion in columen "Encode_main_type". This will add two new columens "Encode_main_cluster" and "Encode_main_cluster_revised".

usage: scRNA_add_cell_type_immune_icell8.py [-h] [--meta META_CLUSTER]
                                            [--black_out BLACK_OUT]
                                            [--replace REPLACE] [-v]

optional arguments:
  -h, --help            show this help message and exit
  --meta META_CLUSTER
  --black_out BLACK_OUT
  --replace REPLACE
  -v

Input:
Cell.ID orig.ident      nCount_RNA      nFeature_RNA
FEL001_C25287_104812_C36_R47    FEL001  12137   845
FEL001_C25287_104812_C49_R49    FEL001  83497   1074
FEL002_C25288_104770_C14_R05    FEL002  18245   784
FEL002_C25288_104770_C22_R07    FEL002  19053   935


Output:
A cross table comparing "Encode_main_cluster" and "Encode_main_cluster_revised":
FEL001010_icell8.immune.metadata.revised_anno.crosstab.txt

A revised annotation file with two new columen "Encode_main_cluster" and "Encode_main_cluster_revised":
FEL001010_icell8.immune.metadata.revised_anno.txt


    '''
    print message


#Cell.ID orig.ident      nCount_RNA	nFeature_RNA Encode_main_cluster     seurat_clusters 
#FEL013_FEL013_E_TTCACCGAGGTATCTC        FEL013  2267    997 Encode_main_cluster     seurat_clusters
def read_meta_with_cluster(infile, cluster, cluster_anno, black_out, epithelial):
    clusters     = defaultdict(lambda : str())
    sample_anno  = defaultdict(lambda : defaultdict(lambda : int()))
    linen     = 0
    output = re.sub(r'.txt', r'.revised_anno.txt', infile)
    ofile  = open(output, 'w')
    cluster_index = 0
    anno_index    = 0
    sample_index  = 0
    encode_index  = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            if len(line) > 2 and line.startswith(r'Cell.ID'):
                for index, value in enumerate(unit):
                    if value == cluster:
                        cluster_index = index
                    elif value == "Encode_main_cluster_revised":
                        anno_index    = index
                    elif value == "Sample":
                        sample_index  = index
                    elif value == 'Encode_main_type':
                        encode_index  = index
                print >> ofile, '\t'.join(unit)
            elif len(line) > 2:
                #print line
                #print unit[0], cluster_index, unit[cluster_index], anno_index, cluster_anno[unit[cluster_index]]
                if black_out.has_key(unit[0]):
                    #print 'Unclassified'
                    unit[anno_index] = 'Unclassified'
                    print >> ofile, '\t'.join(unit)
                else:
                    #print cluster_anno[unit[cluster_index]]
                    if int(epithelial) == 1:
                        if unit[encode_index] in ["B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Endothelial cells", "Eosinophils", "Macrophages", "Monocytes", "Neutrophils", "NK cells"]:
                            unit[anno_index] = 'Low-quality cells'
                        else:
                            unit[anno_index] = cluster_anno[unit[cluster_index]]
                        print >> ofile, '\t'.join(unit)
                    else:
                        unit[anno_index] = cluster_anno[unit[cluster_index]]
                        print >> ofile, '\t'.join(unit)
    ofile.close()

def read_meta_with_cluster_pd(infile, black_out):
    clusters     = defaultdict(lambda : str())
    sample_anno  = defaultdict(lambda : defaultdict(lambda : int()))
    output = re.sub(r'.txt', r'.revised_anno.txt', infile)
    output1 = re.sub(r'.txt', r'.revised_anno.crosstab1.txt', infile)
    output1x = re.sub(r'.txt', r'.revised_anno.crosstab1.xlsx', infile)
    output2 = re.sub(r'.txt', r'.revised_anno.crosstab2.txt', infile)
    output2x = re.sub(r'.txt', r'.revised_anno.crosstab2.xlsx', infile)

    file_df = pd.read_table(infile, header=0)
    col_new_cell_type1 = []
    col_new_cell_type2 = []
    for i in range(file_df.shape[0]):
        if black_out.has_key(str(file_df['Cell.ID'][i])):
            col_new_cell_type1.append('Unclassified')
            col_new_cell_type2.append('Unclassified')
        elif file_df['Encode_main_type'][i] in ["B-cells", "CD4+ T-cells", "CD8+ T-cells", "HSC", "DC", "Macrophages", "Monocytes"]:
            if file_df['Encode_main_type'][i] in ["B-cells"]:
                col_new_cell_type1.append('B cells')
                col_new_cell_type2.append('B cells') 
            elif file_df['Encode_main_type'][i] in ["CD4+ T-cells", "CD8+ T-cells"]:
                col_new_cell_type1.append('T cells')
                col_new_cell_type2.append(re.sub(r'T-cells', r'T cells', file_df['Encode_main_type'][i])) 
            elif file_df['Encode_main_type'][i] in ["DC", "Macrophages", "Monocytes"]:
                col_new_cell_type1.append('Macrophages')
                col_new_cell_type2.append(file_df['Encode_main_type'][i])
            else:
                col_new_cell_type1.append(file_df['Encode_main_type'][i])
                col_new_cell_type2.append(file_df['Encode_main_type'][i])
        else:
            col_new_cell_type1.append('Low-quality cells')
            col_new_cell_type2.append('Low-quality cells')
        #print col_new_cell_type1
    file_df['Encode_main_cluster'] = col_new_cell_type1
    file_df['Encode_main_cluster_revised'] = col_new_cell_type2
    file_df.to_csv(output, sep="\t", index=False)

    com = pd.crosstab(file_df['Encode_main_cluster'], file_df['Encode_main_cluster_revised'])
    com.to_csv(output1, sep="\t", index=True)
    com.to_excel(output1x, index=True)
    com = pd.crosstab(file_df['Encode_main_cluster'], file_df['Encode_main_type'])
    com.to_csv(output2, sep="\t", index=True)
    com.to_excel(output2x, index=True)


def read_black_out(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Cell.ID'):
                unit = re.split(r'\t',line)
                data[unit[0]] = 1
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta')
    parser.add_argument('--black_out')
    parser.add_argument('--replace')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.meta) > 0
    except:
        usage()
        sys.exit(2)

    if not args.replace:
        args.replace = 0

    black_out    = defaultdict(lambda : str())
    if args.black_out:
        black_out = read_black_out(args.black_out)
    #cluster_anno = read_cluster_anno(args.cluster_anno)
    if int(args.replace) == 0:
        read_meta_with_cluster_pd(args.meta, black_out) 
    else:
        read_meta_with_cluster(args.meta, black_out)

if __name__ == '__main__':
    main()

