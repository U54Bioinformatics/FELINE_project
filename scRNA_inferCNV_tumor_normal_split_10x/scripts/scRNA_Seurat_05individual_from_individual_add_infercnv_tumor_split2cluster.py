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
python scRNA_Seurat_05individual_from_individual_add_infercnv_tumor_split2cluster.py --infercnv_predict ../InferCNV_heatmap/FEL026P140_HMM_infer.clust2.anno.txt --meta_cluster FEL026P140_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt --tumor_cluster Cluster2

Use infercnv (two methods: split into two cluster or tirosh's correlation methods) to revise cell type anntoations.
1. annotate cancer cells: epitheliel cells that have CNA
2. remove normal cells that have CNA and assigned to unclassified. These cells might be doublets or have other problems in annotations. Leave these right now.

--infercnv_predict: infercnv results of split into two clusters
--infercnv_tirosh:  infercnv results of tirosh's correlation methods
--meta_cluster: annotation file
--tumor_cluster: clusters that are assigned as tumor based on heatmapcomplex results

    '''
    print message


#Cell.ID orig.ident      nCount_RNA	nFeature_RNA Encode_main_cluster     seurat_clusters 
#FEL013_FEL013_E_TTCACCGAGGTATCTC        FEL013  2267    997 Encode_main_cluster     seurat_clusters
def read_meta_with_cluster(infile, cell_infercnv, cell_tirosh):
    clusters     = defaultdict(lambda : str())
    sample_anno  = defaultdict(lambda : defaultdict(lambda : int()))
    linen     = 0
    output = re.sub(r'.txt', r'.infercnv_tumor_normal.txt', infile)
    file_df = pd.read_table(infile, header=0)
    col_new_cell_type = []
    col_infercnv_split = []
    col_CNA = []
    for i in range(file_df.shape[0]):
        #print file_df[0][i]
        if cell_infercnv.has_key(file_df['Cell.ID'][i]) or cell_tirosh.has_key(file_df['Cell.ID'][i]):
            if file_df['Encode_main_cluster'][i] == 'Epithelial cells':
                col_new_cell_type.append('Cancer cells')
                col_infercnv_split.append('Tumor')
            else:
                col_new_cell_type.append('Unclassified')
                col_infercnv_split.append('Unclassified')
            col_CNA.append('CNA')
        else:
            col_new_cell_type.append(file_df['Encode_main_cluster'][i])
            col_infercnv_split.append('Normal')
            col_CNA.append('noCNA')

    file_df['Infercnv_CNA'] = col_CNA
    file_df['Infercnv_split'] = col_infercnv_split
    file_df['Encode_main_cluster_revised'] = col_new_cell_type
    file_df.to_csv(output, sep="\t", index=False) 


#Cell.ID Sample  nCount_RNA      Total_Reads_k   nFeature_RNA    S       G2M     Phase   percent.mt      Cell_type       seurat_clusters HMM_inferk2
#FEL026_M_AAACCCAGTGCTTCAA       FEL026P140_M    10275   10.275  4117    -0.076663097263 -0.0235519048539        G1      0.5     Epithelial cells        2       Cluster2
def read_tumor_normal(infile, clusters):
    data = defaultdict(lambda : int())
    tumor_clusters = re.split(r',', clusters)
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        #print file_df['HMM_inferk2'][i]
        cell = file_df['Cell.ID'][i]
        unit = re.split(r'_', file_df['Cell.ID'][i])
        if len(unit) >= 4:
            cell = '_'.join(unit[1:])
        else:
            cell = file_df['Cell.ID'][i]

        if file_df['HMM_infer'][i] in tumor_clusters:
            data[cell] = 1
    return data

#Cell.ID CNA.signal      CNA.corr
#FEL026_M_TTTGATCTCTCCTACG       1.65185886295928        0.930368313857358
#FEL026_M_TTTCATGGTCCACTTC       1.65185886295928        0.930368313857358
#FEL026_M_TTGTTTGCACGTAGTT       1.65185886295928        0.930368313857358
def read_cancer_cells_tirosh(infile):
    data = defaultdict(lambda : int())
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        #print file_df['HMM_inferk2'][i]
        data[file_df['Cell.ID'][i]] = 1
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input')
    parser.add_argument('--infercnv_tirosh')
    parser.add_argument('--infercnv_predict')
    parser.add_argument('--meta_cluster')
    parser.add_argument('--tumor_cluster')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.meta_cluster) > 0
    except:
        usage()
        sys.exit(2)


    #cell_infercnv = read_tumor_normal(args.infercnv_predict)
    cell_infercnv = defaultdict(lambda : int())
    cell_tirosh   = defaultdict(lambda : int())
   
    if args.infercnv_predict:
        cell_infercnv = read_tumor_normal(args.infercnv_predict, args.tumor_cluster)
    if args.infercnv_tirosh:
        cell_tirosh   = read_cancer_cells_tirosh(args.infercnv_tirosh)
    read_meta_with_cluster(args.meta_cluster, cell_infercnv, cell_tirosh) 

if __name__ == '__main__':
    main()

