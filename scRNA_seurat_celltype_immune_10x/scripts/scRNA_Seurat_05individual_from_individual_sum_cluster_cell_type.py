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
python Sum_sample_cell_type_anno.py --cluster_anno FEL012016_10x_SingleR_cluster_table_Blueprint.txt --sample_cluster FEL012016_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.sample_matrix.txt


    '''
    print message


#unify cell types as some cells are annoated more specific than others. Hard to compare between samples/patients
def cell_type_sort(in_type):
    cell_type = {
        "CD4+ T-cells"  : "CD4+ T cells",
        "CD4+ T cells"  : "CD4+ T cells",
        "CD4+ Tcm"      : "CD4+ T cells",
        "CD4+ Tem"      : "CD4+ T cells",
        "CD8+ T-cells"  : "CD8+ T cells",
        "CD8+ T cells"  : "CD8+ T cells",
        "CD8+ Tcm"      : "CD8+ T cells",
        "CD8+ Tem"      : "CD8+ T cells",
        "naive B-cells" : "B cells",
        "Class-switched memory B-cells" : "B cells"
        #"Chondrocytes"  : "Fibroblasts",
        #"Smooth muscle" : "Fibroblasts"
    }
    out_type = in_type
    if cell_type.has_key(in_type):
        out_type = cell_type[in_type]
        
    return out_type

#FEL012016_10x_SingleR_cluster_table_Blueprint.txt
#0       1       2       3       4       5       6       7       8       9       10      11      12      13      14      15      16      17      18
#Adipocytes      2       0       0       6       13      1       20      1       2       12      2       0       9       0       21      7       0       0       46
#Astrocytes      2       0       5       0       2       1       4       22      0       29      3       0       0       0       0       0       0       9       0
#CD4+ T-cells    9       6       0       21      5       0       2       0       4       1       1       0       3       0       2       15      0       0       0
def read_table_Blueprint(infile):
    
    clusters     = defaultdict(lambda : str())
    cluster_anno = defaultdict(lambda : int())
    cluster_max  = defaultdict(lambda : int())
    linen     = 0
    output = re.sub(r'.txt', r'.cluster_anno.txt', infile)
    ofile  = open(output, 'w') 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                if linen == 0:
                    print line
                    for index, cluster in enumerate(unit):
                        cluster = re.sub(r'"', r'', cluster)
                        clusters[str(index)] = cluster
                        print index, cluster
                else:
                    print line
                    cell_type = ''
                    for index, value in enumerate(unit):
                        if int(index) >= 1:
                            cluster_index = int(index)
                            cluster_id    = clusters[str(cluster_index)]
                            if int(value) > int(cluster_max[str(cluster_id)]):
                                cluster_anno[str(cluster_id)] = cell_type
                                cluster_max[str(cluster_id)]  = value
                        elif int(index) == 0:
                            cell_type = cell_type_sort(value)
                linen += 1

    print >> ofile, "Cluster\tMax\tCellType"
    for clust in sorted(cluster_anno.keys(), key=int):
        print >> ofile, '%s\t%s\t%s' %(clust, cluster_max[str(clust)], cluster_anno[str(clust)])
        #print >> ofile, '%s\t%s' %(clust, cluster_anno[str(clust)][0], cluster_anno[str(clust)][1])
    ofile.close()
    #return cluster_anno


#0       1       2       3       4       5       6       7       8       9       10      11      12      13      14
#FEL011_E        0       0       0       0       10      0       0       0       0       6       0       7       1
#FEL011_M        0       0       0       0       2       0       0       0       4       13      0       323     3
def read_sample_matrix(infile, cluster_anno):
    clusters     = defaultdict(lambda : str())
    sample_anno  = defaultdict(lambda : defaultdict(lambda : int()))
    linen     = 0
    output = re.sub(r'.txt', r'.sample_anno.txt', infile)
    ofile  = open(output, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                linen += 1
                if linen == 1:
                    #print line
                    for index, cluster in enumerate(unit):
                        cluster = re.sub(r'"', r'', cluster)
                        clusters[str(index)] = cluster
                        #print index, cluster
                else:
                    sample = ''
                    for index, value in enumerate(unit):
                        #print index, value
                        if int(index) >= 1:
                            cluster_index = int(index) - 1
                            cluster_id    = clusters[str(cluster_index)]
                            #print '%s\t%s' %(sample, cluster_anno[str(cluster_id)])
                            sample_anno[sample][cluster_anno[str(cluster_id)]] += int(value)
                        elif int(index) == 0:
                            sample = value

    print >> ofile, 'Sample\tTotal cells\tEpithelial percent\t%s' %('\t'.join(sorted(set(cluster_anno.values()))))
    for sample in sorted(sample_anno.keys()):
        temp_array = []
        total      = 0
        epithelial = 0
        for cell_type in sorted(sample_anno[sample]):
            print cell_type
            total += sample_anno[sample][cell_type]
            if cell_type == 'Epithelial cells':
                epithelial += sample_anno[sample][cell_type]
            temp_array.append(str(sample_anno[sample][cell_type]))
        epi_cell_rate = int(100*float(epithelial)/float(total))
        print >> ofile, '%s\t%s\t%s%s\t%s' %(sample, str(total), str(epi_cell_rate), "%", '\t'.join(temp_array))
    ofile.close()    

#Cluster Max     CellType
#0       2229    Epithelial cells
#1       1647    Epithelial cells
#2       1375    Epithelial cells
#3       443     Fibroblasts
def read_cluster_anno(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Cluster'):
                unit = re.split(r'\t',line)
                print unit[0], unit[2]
                data[str(unit[0])] = unit[2]
    return data

def compare_table(infile, col1, col2):
    outfile = re.sub(r'.txt', '', infile)
    outfile = '%s.crosstab_cluster_celltype.txt' %(outfile)
    file_df = pd.read_table(infile, header=0)
    com = pd.crosstab(file_df[col1], file_df[col2])
    com.to_csv(outfile, sep="\t", index=True)
    return outfile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta')
    parser.add_argument('--cluster_col')
    parser.add_argument('--celltype_col')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.meta) > 0
    except:
        usage()
        sys.exit(2)


    if not args.cluster_col:
        args.cluster_col = 'RNA_snn_res.0.8'

    if not args.celltype_col:
        args.celltype_col = 'Encode_main_cluster_revised'

    # generate a cross table for cluster and celltype from meta table
    cluster_celltype = compare_table(args.meta, args.celltype_col, args.cluster_col)
    print 'cluster celltype cross table: %s' %(cluster_celltype)
    #
    cluster_anno_file = re.sub(r'.txt', r'.cluster_anno.txt', cluster_celltype)
    if not os.path.exists(cluster_anno_file):
        cluster_anno_file = read_table_Blueprint(cluster_celltype)
    #read cluster annotation from file, so we can revise the annotation through manual inspection
    #cluster_anno      =  read_cluster_anno(cluster_anno_file)
    #read_sample_matrix(args.sample_cluster, cluster_anno)

if __name__ == '__main__':
    main()

