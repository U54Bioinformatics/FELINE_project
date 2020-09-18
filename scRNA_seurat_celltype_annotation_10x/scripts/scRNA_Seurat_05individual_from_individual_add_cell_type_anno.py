#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
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
    sample_index  = 0 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            if len(line) > 2 and line.startswith(r'Cell'):
                for index, value in enumerate(unit):
                    if value == "integrated_snn_res.0.8":
                        cluster_index = index
                    elif value == "Encode_main_cluster":
                        anno_index    = index
                    elif value == "Sample":
                        sample_index = index
                unit.append('Encode_main_cluster_revised')
                print >> ofile, '\t'.join(unit)
            elif len(line) > 2:
                #print line
                #print cluster_index, unit[cluster_index]
                anno_celltype = cluster_anno[unit[cluster_index]]
                unit.append(anno_celltype)
                print >> ofile, '\t'.join(unit)
    ofile.close()


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--cluster_anno')
    parser.add_argument('--meta_cluster')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.cluster_anno) > 0 and len(args.meta_cluster) > 0
    except:
        usage()
        sys.exit(2)


    cluster_anno = read_cluster_anno(args.cluster_anno)
    read_meta_with_cluster(args.meta_cluster, cluster_anno) 

if __name__ == '__main__':
    main()

