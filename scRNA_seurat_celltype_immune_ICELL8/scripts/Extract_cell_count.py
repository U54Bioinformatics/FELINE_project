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
python Extract_cell_count.py --cell test.cells.id.txt --count FEL011024_10x_gene_symbols.filtered.counts.txt


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

#Gene.ID FEL011P138_FEL011_M_AAACGAACACAAGTGG    FEL011P138_FEL011_M_AAACGCTGTTAAGAAC    FEL011P138_FEL011_M_AAAGGTATCAGCCTCT    FEL011P138_FEL011_M_AAAGTCCCAGCGAGTA    FEL011P138_FE
#LINC00115       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       1       0       0       0       0
#LINC02593       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
def read_count(infile, prefix, cells):
    cell_index = defaultdict(lambda : int())
    file_df = pd.read_table(infile, header=0)
    cell_list = ["Gene.ID"]
    cell_list.extend(sorted(cells.keys()))
    file_df_subset = file_df[cell_list]
    outfile = '%s.subset.txt' %(infile)
    file_df_subset.to_csv(outfile, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--count')
    parser.add_argument('--cell')
    parser.add_argument('--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    try:
        len(args.count) > 0 and len(args.cell) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = 'Sample'

    cells = read_cell_id(args.cell)
    read_count(args.count, args.project, cells)

if __name__ == '__main__':
    main()

