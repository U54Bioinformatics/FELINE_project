#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import re
import os
import argparse
import loompy

def usage():
    test="name"
    message='''
/home/jichen/software/BETSY/install/envs/scRNA_velocyto/bin/python scRNA_RNA_velocity_merge_loom_files.py --bamdir ../Preprocess/VG_CAMA1_D11_ALL_count03/ --sample samples.list

    '''
    print(message)

def read_sample_file(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=None)
    file_df.columns = ["Sample"]
    for i in range(file_df.shape[0]):
        #print file_df[0][i]
        data[file_df['Sample'][i]] = 1
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bamdir')
    parser.add_argument('--sample')
    parser.add_argument('--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.bamdir) > 0 and len(args.sample) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = 'VG_CAMA1_D11_ALL'

    samples = read_sample_file(args.sample)
    #print(samples)
    loomfiles = []
    for sample in sorted(samples.keys()):
        loomfile = '%s/%s/velocyto/%s.loom' %(args.bamdir, sample, sample)
        print(loomfile)
        ds = loompy.connect(loomfile)
        print(ds.shape)
        loomfiles.append(loomfile)

    loompy.combine(loomfiles, '%s.loom' %(args.output), key="Accession")
    ds = loompy.connect('%s.loom' %(args.output))
    print("Merge loom files")
    print (ds.shape)

if __name__ == '__main__':
    main()

