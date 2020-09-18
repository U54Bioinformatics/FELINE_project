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

Reformat icell8 cell type annotation. This resulting has a broader cell type annoation which need to be updated and finalized. The previous results has missing information.

    '''
    print message


#Cell.ID orig.ident      nCount_RNA	nFeature_RNA Encode_main_cluster     seurat_clusters 
#FEL013_FEL013_E_TTCACCGAGGTATCTC        FEL013  2267    997 Encode_main_cluster     seurat_clusters
def read_meta_with_cluster(infile, cell_anno):
    output = re.sub(r'.txt', r'.revised_anno.txt', infile)
    output1 = re.sub(r'.txt', r'.revised_anno.crosstab1.txt', infile)
    output1x = re.sub(r'.txt', r'.revised_anno.crosstab1.xlsx', infile)
    output2 = re.sub(r'.txt', r'.revised_anno.crosstab2.txt', infile)
    output2x = re.sub(r'.txt', r'.revised_anno.crosstab2.xlsx', infile)


    file_df = pd.read_table(infile, header=0)
    header  = ["Cell.ID", "orig.ident", "nCount_RNA", "nFeature_RNA", "Sample", "Cell.orig", "percent.mt", "S.Score", "G2M.Score", "Phase", "hpca_main_type", "Encode_main_type"]
    file_df = file_df[header]
    col_new_cell_type1 = []
    col_new_cell_type2 = []
    for i in range(file_df.shape[0]):
        if cell_anno.has_key(str(file_df['Cell.ID'][i])):
            col_new_cell_type1.append(cell_anno[file_df['Cell.ID'][i]][0])
            col_new_cell_type2.append(cell_anno[file_df['Cell.ID'][i]][1])
        else:
            col_new_cell_type1.append('NA')
            col_new_cell_type2.append('NA')

    file_df['seurat_clusters'] = col_new_cell_type1
    file_df['Encode_main_cluster'] = col_new_cell_type2
    file_df.to_csv(output, sep="\t", index=False)

    com = pd.crosstab(file_df['Encode_main_cluster'], file_df['seurat_clusters'])
    com.to_csv(output1, sep="\t", index=True)
    com.to_excel(output1x, index=True)
    com = pd.crosstab(file_df['Encode_main_cluster'], file_df['Encode_main_type'])
    com.to_csv(output2, sep="\t", index=True)
    com.to_excel(output2x, index=True)
 

#Cell.ID RNA_snn_res.0.8 Encode_main_cluster_revised
def read_cell_anno(infile):
    data = defaultdict(lambda : int())
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        data[file_df['Celltype'][i]] = 1
    
    for ct in sorted(data.keys()):
        ct_name = re.sub(r' ', r'_', ct)
        outfile = re.sub(r'.txt', r'.%s.txt' %(ct_name), infile)
        file_df_ct = file_df[file_df['Celltype'] == ct]
        file_df_ct.to_csv(outfile, sep="\t", index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.meta) > 0
    except:
        usage()
        sys.exit(2)


    read_cell_anno(args.meta)

if __name__ == '__main__':
    main()

