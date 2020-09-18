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

#../Infercnv/FEL011P138_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt > FEL011P138_cell_metadata.txt
#Cell.ID Sample  Total.Reads     Total.Reads..Thousands. Expressed.Features      S       G2M     Phase   Percent.Mitochondria    Encode_main_cluster     seurat_clusters
#Cell.ID orig.ident      nCount_RNA      nFeature_RNA    Sample  Cell.orig       percent.mt      S.Score G2M.Score       Phase   old.ident       RNA_snn_res.0.5 hpca_main_type  Encode_main_type     hpca_type       Encode_type     Encode_main_cluster     seurat_clusters
def read_anno_file(infile, output):
    patient = re.split('_', os.path.split(infile)[1])[0]
    outfile = '%s_cell_metadata.txt' %(patient)
    if output:
        outfile = output
    #header  = ["Cell.ID", "Sample", "Total.Reads", "Total.Reads..Thousands.", "Expressed.Features", "S", "G2M", "Phase", "Percent.Mitochondria", "Encode_main_cluster", "seurat_clusters"]
    header  = ["Cell.ID", "Sample", "nCount_RNA", "nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "Phase", "percent.mt", "Encode_main_cluster", "seurat_clusters"]
    file_df = pd.read_table(infile, header=0)
    file_df_sub = file_df[header]
    file_df_sub.to_csv(outfile, sep="\t", index=False)

#../Seurat/FEL028P131_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.infercnv_tumor_normal.txt > FEL028P131_cell_metadata_subclone.txt
#Cell.ID Sample  Total.Reads     Total.Reads..Thousands. Expressed.Features      S       G2M     Phase   Percent.Mitochondria    Encode_main_cluster     seurat_clusters
def read_anno_file_4_subclone(infile, output):
    patient = re.split('_', os.path.split(infile)[1])[0]
    outfile = '%s_cell_metadata_subclone.txt' %(patient)
    if output:
        outfile = output
    header  = ["Cell.ID", "Sample", "Total.Reads", "Total.Reads..Thousands.", "Expressed.Features", "S", "G2M", "Phase", "Percent.Mitochondria", "Encode_main_cluster", "seurat_clusters", "Infercnv_split"]
    file_df = pd.read_table(infile, header=0)
    file_df_sub = file_df[header]
    file_df_sub = file_df_sub[file_df_sub["Infercnv_split"] == 'Tumor']
    file_df_sub.to_csv(outfile, sep="\t", index=False)
 

def main():
   
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta')
    parser.add_argument('--output')
    parser.add_argument('--subclone', action='store_true')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.meta) > 0
    except:
        usage()
        sys.exit(2)

    if args.subclone:
        read_anno_file_4_subclone(args.meta, args.output)
    else:
        read_anno_file(args.meta, args.output)

if __name__ == '__main__':
    main()
