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
python Add_cell_type_to_meta.py --cluster_anno FEL012016_10x_SingleR_cluster_table_Blueprint.cluster_anno.txt --meta_cluster FEL012016_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.txt


    '''
    print message


#Cell.ID orig.ident      nCount_RNA	nFeature_RNA Encode_main_cluster     seurat_clusters 
#FEL013_FEL013_E_TTCACCGAGGTATCTC        FEL013  2267    997 Encode_main_cluster     seurat_clusters
def read_meta_with_cluster(infile, cluster, cluster_anno, black_out, epithelial, replace_anno):
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
                    elif value == replace_anno:
                        anno_index    = index
                    elif value == "Sample":
                        sample_index  = index
                    #elif value == "Encode_main_cluster"
                    elif value == 'Celltype':
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
                        if cluster_anno[unit[cluster_index]] in ['Macrophages', 'Monocytes', "DC"]:
                            if unit[encode_index] in ['Macrophages', 'Monocytes', "DC"]:
                                unit[anno_index] = cluster_anno[unit[cluster_index]]
                            else:
                                unit[anno_index] = 'Low-quality cells'
                        elif cluster_anno[unit[cluster_index]] in ['T cells', 'CD4+ T cells', 'CD8+ T cells', 'NK cells', "Tregs"]:
                            if unit[encode_index] in ['T cells', 'CD4+ T cells', 'CD8+ T cells', 'NK cells', 'Tregs']:
                                unit[anno_index] = cluster_anno[unit[cluster_index]]
                            else: 
                                unit[anno_index] = 'Low-quality cells'
                        elif cluster_anno[unit[cluster_index]] in ['B cells', 'Plasma cells']:
                            #print cluster_anno[unit[cluster_index]]
                            if unit[encode_index] in ['B-cells', 'Plasma-cells', 'B cells', 'Plasma cells']:
                                unit[anno_index] = cluster_anno[unit[cluster_index]]
                            else:
                                unit[anno_index] = 'Low-quality cells'
                        elif cluster_anno[unit[cluster_index]] in ['Endothelial_cells', 'Endothelial cells', 'Vas-Endo', 'Lym-Endo']:
                            if unit[encode_index] in ['Endothelial_cells', 'Endothelial cells', 'Vas-Endo', 'Lym-Endo']:
                                unit[anno_index] = cluster_anno[unit[cluster_index]]
                            else:
                                unit[anno_index] = 'Low-quality cells'
                        else:
                            unit[anno_index] = 'Low-quality cells'
                        #unit[anno_index] = cluster_anno[unit[cluster_index]]
                        print >> ofile, '\t'.join(unit)
    ofile.close()

def read_meta_with_cluster_pd(infile, cluster, cluster_anno, black_out, epithelial):
    clusters     = defaultdict(lambda : str())
    sample_anno  = defaultdict(lambda : defaultdict(lambda : int()))
    output = re.sub(r'.txt', r'.revised_anno.txt', infile)
    file_df = pd.read_table(infile, header=0)
    col_new_cell_type1 = []
    for i in range(file_df.shape[0]):
        #print file_df[0][i]
        #print file_df[cluster][i]
        if cluster_anno.has_key(str(file_df[cluster][i])):
            if epithelial == 1:
                if file_df['Encode_main_type'][i] in ["B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Endothelial cells", "Eosinophils", "Macrophages", "Monocytes", "Neutrophils", "NK cells"]:
                    col_new_cell_type1.append('Low-quality cells')
                else:
                    col_new_cell_type1.append(cluster_anno[str(file_df[cluster][i])])
            else:
                col_new_cell_type1.append(cluster_anno[str(file_df[cluster][i])])
        elif black_out.has_key(str(file_df[cluster][i])):
            col_new_cell_type1.append('Unclassified')
        else:
            col_new_cell_type1.append('NA')
        #print col_new_cell_type1
    file_df['Encode_main_cluster_revised'] = col_new_cell_type1
    file_df.to_csv(output, sep="\t", index=False)

#Cluster CellType
#0       Epithelial cells
#1       Fibroblasts
#2       Fibroblasts
def read_cluster_anno(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        data[str(file_df['Cluster'][i])] = file_df['CellType'][i]
    return data

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
    parser.add_argument('--cluster')
    parser.add_argument('--cluster_anno')
    parser.add_argument('--replace_anno')
    parser.add_argument('--meta_cluster')
    parser.add_argument('--black_out')
    parser.add_argument('--epithelial')
    parser.add_argument('--replace')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.cluster_anno) > 0 and len(args.meta_cluster) > 0
    except:
        usage()
        sys.exit(2)

    if not args.cluster:
        args.cluster = 'SCT_snn_res.0.8'

    if not args.replace:
        args.replace = 0

    if not args.epithelial:
        args.epithelial = 0

    if not args.replace_anno:
        args.replace_anno = 'Encode_main_cluster_revised'

    black_out    = defaultdict(lambda : str())
    if args.black_out:
        black_out = read_black_out(args.black_out)
    cluster_anno = read_cluster_anno(args.cluster_anno)
    if int(args.replace) == 0:
        read_meta_with_cluster_pd(args.meta_cluster, args.cluster, cluster_anno, black_out, args.epithelial) 
    else:
        read_meta_with_cluster(args.meta_cluster, args.cluster, cluster_anno, black_out, args.epithelial, args.replace_anno)

if __name__ == '__main__':
    main()

