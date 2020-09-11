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
    file_obs = glob.glob("./%s/*.observations.txt" %(indir))
    file_png = glob.glob("./%s/*.png" %(indir))
    file_in  = file_obs + file_png
    #file_in  = file_png
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
        args.prefix = "PREFIX"


    filenames = get_file_names(args.infercnv_outdir)
    for fn in sorted(filenames):
        file_old = '%s/%s' %(args.infercnv_outdir, fn)
        file_new = '%s/%s_%s' %(args.infercnv_outdir, args.prefix, fn)
        print file_old
        print file_new
        os.system('cp %s %s' %(file_old, file_new))

if __name__ == '__main__':
    main()

