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


def convert_patient_id(id_in):
    patients = {
         "FEL001" : "FEL001P120",
         "FEL002" : "FEL002P106",
         "FEL003" : "FEL003P112",
         "FEL004" : "FEL004P114",
         "FEL005" : "FEL005P133",
         "FEL006" : "FEL006P111",
         "FEL007" : "FEL007P127",
         "FEL008" : "FEL008P108",
         "FEL009" : "FEL009P136",
         "FEL010" : "FEL010P104"
    }
    id_out = patients[id_in]
    return id_out

#Cell.ID orig.ident      nCount_RNA	nFeature_RNA Encode_main_cluster     seurat_clusters 
#FEL013_FEL013_E_TTCACCGAGGTATCTC        FEL013  2267    997 Encode_main_cluster     seurat_clusters
def finalize_celltype_anno(infile, infercnv_anno, immune_anno):
    output = re.sub(r'.txt', r'.finalized.txt', infile)
    output_short = re.sub(r'.txt', r'.finalized_short.txt', infile)
    output1 = re.sub(r'.txt', r'.finalized.crosstab1.txt', infile)
    output1x = re.sub(r'.txt', r'.finalized.crosstab1.xlsx', infile)
    output2 = re.sub(r'.txt', r'.finalized.crosstab2.txt', infile)
    output2x = re.sub(r'.txt', r'.finalized.crosstab2.xlsx', infile)
    output_sample_anno = re.sub(r'.txt', r'.finalized.sample_celltype.txt', infile)
    output_sample_annox = re.sub(r'.txt', r'.finalized.sample_celltype.xlsx', infile)
    output_sample_anno1 = re.sub(r'.txt', r'.finalized.sample_celltype_subtype.txt', infile)
    output_sample_anno1x = re.sub(r'.txt', r'.finalized.sample_celltype_subtype.xlsx', infile)

    file_df = pd.read_table(infile, header=0)
    col_new_cell_type1 = []
    col_new_cell_type2 = []
    col_new_cell_type3 = []
    col_new_cell_type4 = []
    for i in range(file_df.shape[0]):
        # Infercnv_CNA
        cna = 0
        if infercnv_anno.has_key(str(file_df['Cell.ID'][i])):
            col_new_cell_type3.append(infercnv_anno[file_df['Cell.ID'][i]][0])
            if infercnv_anno[file_df['Cell.ID'][i]][0] == 'CNA':
                cna = 1
        else:
            col_new_cell_type3.append('NA')
       
        # platform
        col_new_cell_type4.append("ICELL8")
  
        # Sample
        file_df['Sample'].iloc[i,] = re.sub(r'%s' %(file_df['orig.ident'][i]), r'%s' %(convert_patient_id(file_df['orig.ident'][i])), file_df['Sample'][i])

        # Celltype and Celltype_subtype
        if file_df['Encode_main_cluster'][i] in ["Epithelial cells"]:
            if cna == 0:
                col_new_cell_type1.append("Normal epithelial cells")
                col_new_cell_type2.append("Normal epithelial cells")
            else:
                col_new_cell_type1.append("Cancer cells")
                col_new_cell_type2.append("Cancer cells")
        elif file_df['Encode_main_cluster'][i] in ["Immune"]:
            if cna == 0:
                col_new_cell_type1.append(immune_anno[file_df['Cell.ID'][i]][0])
                col_new_cell_type2.append(immune_anno[file_df['Cell.ID'][i]][1])
            else:
                #col_new_cell_type1.append(immune_anno[file_df['Cell.ID'][i]][0])
                col_new_cell_type1.append("Low-quality cells")
                col_new_cell_type2.append("Low-quality cells")
        else:
            if cna == 0:
                col_new_cell_type1.append(file_df['Encode_main_cluster'][i])
                col_new_cell_type2.append(file_df['Encode_main_cluster'][i])
            else:
                #col_new_cell_type1.append(file_df['Encode_main_cluster'][i])
                col_new_cell_type1.append("Low-quality cells")
                col_new_cell_type2.append("Low-quality cells")


    # final
    file_df['Celltype'] = col_new_cell_type1
    file_df['Celltype_subtype'] = col_new_cell_type2
    file_df['Infercnv_CNA'] = col_new_cell_type3
    file_df['Platform'] = col_new_cell_type4
    file_df.to_csv(output, sep="\t", index=False)

    # short
    header  = ["Cell.ID", "orig.ident", "nCount_RNA", "nFeature_RNA", "Sample", "Celltype", "Celltype_subtype", "Infercnv_CNA", "Platform"]
    file_df_short = file_df[header]
    file_df_short.to_csv(output_short, sep="\t", index=False)

    # comparison and summary
    com = pd.crosstab(file_df['Encode_main_cluster'], file_df['Celltype'])
    com.to_csv(output1, sep="\t", index=True)
    com.to_excel(output1x, index=True)
    com = pd.crosstab(file_df['Encode_main_cluster'], file_df['Celltype_subtype'])
    com.to_csv(output2, sep="\t", index=True)
    com.to_excel(output2x, index=True)
    com = pd.crosstab(file_df['Sample'], file_df['Celltype'])
    com.to_csv(output_sample_anno, sep="\t", index=True)
    com.to_excel(output_sample_annox, index=True)
    com = pd.crosstab(file_df['Sample'], file_df['Celltype_subtype'])
    com.to_csv(output_sample_anno1, sep="\t", index=True)
    com.to_excel(output_sample_anno1x, index=True)

   

#Infercnv_CNA    Infercnv_split  Encode_main_cluster_revised
def read_infercnv_anno(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        data[file_df['Cell.ID'][i]] = [file_df['Infercnv_CNA'][i]]
    return data

#Encode_main_cluster     Encode_main_cluster_revised
def read_immune_anno(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        data[file_df['Cell.ID'][i]] = [file_df['Encode_main_cluster'][i], file_df['Encode_main_cluster_revised'][i]]
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta')
    parser.add_argument('--infercnv_anno')
    parser.add_argument('--immune_anno')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.meta) > 0 and len(args.infercnv_anno) > 0 and len(args.immune_anno) > 0
    except:
        usage()
        sys.exit(2)


    infercnv_anno = read_infercnv_anno(args.infercnv_anno)
    immune_anno = read_immune_anno(args.immune_anno)
    finalize_celltype_anno(args.meta, infercnv_anno, immune_anno) 

if __name__ == '__main__':
    main()

