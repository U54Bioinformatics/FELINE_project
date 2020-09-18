#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import re
import os
import argparse

def usage():
    test="name"
    message='''
python Add_cell_type_to_meta.py --cluster_anno FEL012016_10x_SingleR_cluster_table_Blueprint.cluster_anno.txt --meta_cluster FEL012016_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.txt


    '''
    print message


#Cell.ID orig.ident      nCount_RNA	nFeature_RNA Encode_main_cluster     seurat_clusters 
#FEL013_FEL013_E_TTCACCGAGGTATCTC        FEL013  2267    997 Encode_main_cluster     seurat_clusters
def read_meta_with_cluster(infile, cluster_anno):
    clusters     = defaultdict(lambda : str())
    sample_anno  = defaultdict(lambda : defaultdict(lambda : int()))
    linen     = 0
    output = re.sub(r'.txt', r'.revised_anno.txt', infile)
    ofile  = open(output, 'w')
    cluster_index = 0
    anno_index    = 0
    anno_sub_index= 0
    sample_index  = 0
    infercnv_index= 0 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            if len(line) > 2 and line.startswith(r'Cell'):
                for index, value in enumerate(unit):
                    if value == "seurat_clusters":
                        cluster_index = index
                    elif value == "orig.ident":
                        sample_index = index
                    elif value == 'Celltype1':
                        anno_index = index
                    elif value == 'Celltype2':
                        anno_sub_index = index
                    elif value == 'Infercnv_CNA':
                        infercnv_index = index
                unit.append('Celltype')
                unit.append('Celltype_subtype')
                print >> ofile, '\t'.join(unit)
            elif len(line) > 2:
                #print line
                #print cluster_index, unit[cluster_index]
                anno_celltype = ''
                anno_subtype  = ''
                # if patients are new (#44, #45, #46) use new annotation, otherwise use previous annotation
                if unit[sample_index] in ['FEL044', 'FEL045', 'FEL046']:
                    anno_celltype = cluster_anno[unit[cluster_index]]
                    anno_subtype  = anno_celltype
                else:
                    anno_celltype = unit[anno_index]
                    anno_subtype  = unit[anno_sub_index]
                # for epithelial cells, set these with CNA as cancer cells and these wo CNA as normal
                if anno_celltype in ['Normal epithelial cells', 'Cancer cells']:
                    if unit[infercnv_index] == 'CNA':
                        anno_celltype = 'Cancer cells'
                        anno_subtype  = 'Cancer cells'
                    elif unit[infercnv_index] == 'noCNA':
                        anno_celltype = 'Normal epithelial cells'
                        anno_subtype  = 'Normal epithelial cells'
                # for other cells, set these with CNA as Low-quality cells and these wo CNA as annoated
                else:
                    if unit[infercnv_index] == 'CNA':
                        anno_celltype = 'Low-quality cells'
                        anno_subtype  = 'Low-quality cells'
                unit.append(anno_celltype)
                unit.append(anno_subtype)
                print >> ofile, '\t'.join(unit)
    ofile.close()
    return output

#Cluster Max     CellType
#0       1021    Epithelial cells
#1       982     Epithelial cells
#2       815     Epithelial cells
#3       667     Epithelial cells
def read_cluster_anno(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Cluster'):
                unit = re.split(r'\t',line)
                data[unit[0]] = unit[2]
    return data


#extract these column: Cell.ID orig.ident      nCount_RNA      nFeature_RNA    Sample  Celltype        Celltype_subtype        Infercnv_CNA
#Cell.ID orig.ident      nCount_RNA      nFeature_RNA    Sample  Celltype        Celltype_subtype	Infercnv_CNA    Platform
def short_anno_file(infile):
    output = re.sub(r'.txt', r'.short.txt', infile)
    file_df = pd.read_table(infile, header=0)
    header  = ['Cell.ID', 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'Sample', 'Celltype', 'Celltype_subtype', 'Infercnv_CNA', 'Platform']
    file_df = file_df[header]
    file_df.to_csv(output, index=False, sep="\t")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)


    short_anno_file(args.input)

if __name__ == '__main__':
    main()

