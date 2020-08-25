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
python Prepare_celltype_cell_list.py --meta FEL011046_10x_cell_metadata.UMAPcluster.integrated_RPCA.infercnv_anno.scrublet_anno.revised_anno.short.txt --output FEL011046 --platform 10x
 
Generate Cell.ID files for each cell type. The Cell.ID files can be used to extract meta and count from master files
    '''
    print message


#Cell.ID orig.ident      nCount_RNA      nFeature_RNA    Sample  Celltype        Celltype_subtype        Infercnv_CNA
#FEL011_M_AAACGAACACAAGTGG       FEL011  5390    2361    FEL011P138_M    Cancer cells    Cancer cells    CNA
#FEL011_M_AAACGCTGTTAAGAAC       FEL011  5394    2123    FEL011P138_M    Cancer cells    Cancer cells    CNA
def read_cell_anno(infile):
    data = defaultdict(lambda : list())
    file_df = pd.read_table(infile, sep = "\t", header=0)
    for i in range(file_df.shape[0]):
        #print file_df['Encode_main_cluster'][i], file_df['Cell.ID'][i]
        celltype = file_df['Celltype'][i]
        celltype = re.sub(r' ', r'_', celltype)
        data[celltype].append(file_df['Cell.ID'][i])
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta')
    parser.add_argument('--output')
    parser.add_argument('--platform')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.meta) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = 'FEL001010'
 
    if not args.platform:
        args.platform = 'icell8'

    cell_anno = read_cell_anno(args.meta)
    # write cell id of each cell into a file
    for celltype in sorted(cell_anno.keys()):
        celltype = re.sub(r' ', r'_', celltype)
        outfile = '%s_%s_%s.Cell_ID.txt' %(args.output, celltype, args.platform)
        ofile = open(outfile, 'w')
        for cell in sorted(cell_anno[celltype]):
            print >> ofile, cell
            #print cell, celltype
        ofile.close()

if __name__ == '__main__':
    main()

