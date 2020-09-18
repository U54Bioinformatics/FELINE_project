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
python Extract_cell_meta.py --cell test.cells.id.txt --meta FEL011024_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.txt


    '''
    print message


#FEL014P103_FEL014_S_GTAGAAAAGAAGCGGG
def read_cell_id(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            data[unit[0]] = 1
    return data

#Cell.ID orig.ident      nCount_RNA      nFeature_RNA    Sample  Total.Reads     Total.Reads..Thousands. Total.Non.Mitochondria.Reads    Total.Non.Mitochondria.Reads..Thousands.
#FEL022P124_FEL022_S_TGTGTGACACAATGTC    FEL022  3025    1819    FEL022P124_S    3025    3.025   3015    3.015   1819    0.154389457664  0.465553133701  G2M     0.3     FEL022_S_TGTG
#FEL022P124_FEL022_M_AGGAAATAGTGGCCTC    FEL022  4493    2247    FEL022P124_M    4493    4.493   4460    4.46    2247    -0.00673886300505       -0.0552656534381        G1      0.7
def read_meta(infile, prefix, cells):
    cell_index = defaultdict(lambda : int())
    outfile = '%s.subset.txt' %(infile)
    ofile = open(outfile, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            if unit[0] == 'Cell.ID':
                print >> ofile, line
            else:
                if cells.has_key(unit[0]):
                    print >> ofile, line
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta')
    parser.add_argument('--cell')
    parser.add_argument('--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    try:
        len(args.meta) > 0 and len(args.cell) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = 'Sample'

    cells = read_cell_id(args.cell)
    read_meta(args.meta, args.project, cells)

if __name__ == '__main__':
    main()

