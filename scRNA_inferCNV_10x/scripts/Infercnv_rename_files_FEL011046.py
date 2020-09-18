#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob

def usage():
    test="name"
    message='''
python Infercnv_rename_files.py --infercnv_outdir FEL013_10x_counts.all.output_dir_plot --prefix FEL013_10x_all_normal

    '''
    print message


def get_file_names(indir):
    file_in = ["infercnv.preliminary.observations.txt", "infercnv.observations.txt", "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"]
    files = []
    for f in sorted(file_in):
        fn = os.path.split(f)[1]
        files.append(fn)
    return files


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infercnv_outdir')
    parser.add_argument('--prefix')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.infercnv_outdir) > 0
    except:
        usage()
        sys.exit(2)

    if not args.prefix:
        subfolder = os.path.split(os.path.abspath(args.infercnv_outdir))[1]
        args.prefix = re.split(r'_', subfolder)[0]
        #print args.prefix

    filenames = get_file_names(args.infercnv_outdir)
    for fn in sorted(filenames):
        file_old = '%s/%s' %(args.infercnv_outdir, fn)
        file_new = '%s/%s_%s' %(args.infercnv_outdir, args.prefix, fn)
        if os.path.exists(file_old):
            print file_old
            print file_new
            os.system('cp %s %s' %(file_old, file_new))

if __name__ == '__main__':
    main()

