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
python scRNA_Seurat_05individual_from_individual_add_SCINA_anno.py --scina FEL011024_10x_gene_symbols.normalized.counts.SCINA.txt --meta_cluster FEL011024_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt


    '''
    print message


#Cell.ID orig.ident      nCount_RNA	nFeature_RNA Encode_main_cluster     seurat_clusters 
#FEL013_FEL013_E_TTCACCGAGGTATCTC        FEL013  2267    997 Encode_main_cluster     seurat_clusters
def read_meta_with_cluster(infile, cell_anno):
    clusters     = defaultdict(lambda : str())
    sample_anno  = defaultdict(lambda : defaultdict(lambda : int()))
    output = re.sub(r'.txt', r'.scrublet_anno.txt', infile)
    file_df = pd.read_table(infile, header=0)
    col_new_cell_type1 = []
    col_new_cell_type2 = []
    col_new_cell_type3 = []
    for i in range(file_df.shape[0]):
        #print file_df[0][i]
        if cell_anno.has_key(file_df['Cell.ID'][i]):
            col_new_cell_type1.append(cell_anno[file_df['Cell.ID'][i]][0])
            col_new_cell_type2.append(cell_anno[file_df['Cell.ID'][i]][1])
            col_new_cell_type3.append(cell_anno[file_df['Cell.ID'][i]][2])
        else:
            col_new_cell_type1.append('NA')
            col_new_cell_type2.append('NA')
            col_new_cell_type3.append('NA')

    file_df['scrublet_doublet_score'] = col_new_cell_type1
    file_df['scrublet_doublet_call1'] = col_new_cell_type2
    file_df['scrublet_doublet_call2'] = col_new_cell_type3
    file_df.to_csv(output, sep="\t", index=False) 


#Cell.ID scrublet_doublet_score  scrublet_doublet_call1  scrublet_doublet_call2
#FEL011_S_AAACCCAGTGGCGCTT       0.176161262051  False   False
#FEL011_S_AAACCCATCAAACGAA       0.0558035714286 False   False
#FEL011_S_AAACCCATCCTACTGC       0.176161262051  False   False
def read_scrublet_anno(infile):
    data = defaultdict(lambda : list())
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        call1 = 'NA'
        call2 = 'NA'
        #scrublet_doublet_call1
        if file_df['scrublet_doublet_call1'][i]:
            call1 = 'Doublet'
        else:
            call1 = 'Singlet'
        #scrublet_doublet_call2
        if file_df['scrublet_doublet_call2'][i]:
            call2 = 'Doublet'
        else:
            call2 = 'Singlet'
        #score, call1 and call2
        data[file_df['Cell.ID'][i]] = [file_df['scrublet_doublet_score'][i], call1, call2]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input')
    parser.add_argument('--scrublet')
    parser.add_argument('--meta_cluster')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.meta_cluster) > 0
    except:
        usage()
        sys.exit(2)


    cell_anno = read_scrublet_anno(args.scrublet)
    read_meta_with_cluster(args.meta_cluster, cell_anno) 

if __name__ == '__main__':
    main()

