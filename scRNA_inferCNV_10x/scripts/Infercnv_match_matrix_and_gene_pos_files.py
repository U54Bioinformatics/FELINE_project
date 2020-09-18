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
python Match_matrix_and_gene_pos_files.py --project 97OVCZ_10x_counts

    '''
    print message



#DDX11L1 chr1    11869   14412
#WASH7P  chr1    14363   29806
def read_gene_pos(infile):
    data = defaultdict(lambda : str())
    genes= []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                data[unit[0]] = line
                genes.append(unit[0])
    return data, genes


#ENSG00000243485 MIR1302-10
#ENSG00000237613 FAM138A
def read_gene_tsv(infile):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                data.append(unit[1])
    return data


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--project')
    parser.add_argument('-m', '--matrix')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.matrix) > 0
    except:
        usage()
        sys.exit(2)

   
    #genes_in_10x = read_gene_tsv("./filtered_gene_bc_matrices/hg19/genes.tsv")
    genes_in_pos, gene_sorted = read_gene_pos("hg19.RefSeq.NM_pos_unique_sort.txt")

    #ofile = open ("%s.gene_pos.txt" %(args.project), 'w')
    #for g in genes_in_10x:
    #    if genes_in_pos.has_key(g):
    #        print >> ofile, genes_in_pos[g] 
    #ofile.close()

    matrix_raw = args.matrix
    matrix_file= re.sub(r'.raw.matrix', r'.matrix', matrix_raw)
    line_num  = 0
    gene_list = defaultdict(lambda : int())
    ofile = open(matrix_file, 'w')
    with open (matrix_raw, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            line_num += 1
            if line_num == 1:
                print >> ofile, line
            else:
                #print unit[0]
                if genes_in_pos.has_key(unit[0]):
                    print >> ofile, line
                    gene_list[unit[0]] = 1
                else:
                    continue
    ofile.close()

    gene_pos = re.sub(r'.raw.matrix', r'.gene_pos.txt', matrix_raw)
    ofile = open (gene_pos, 'w')
    for g in gene_sorted:
        if gene_list.has_key(g):
            print >> ofile, genes_in_pos[g]
    ofile.close()


if __name__ == '__main__':
    main()
