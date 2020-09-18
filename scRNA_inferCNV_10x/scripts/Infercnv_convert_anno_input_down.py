#!/usr/bin/python
import sys
from collections import defaultdict
import re
import os
import argparse
import glob


def usage():
    test="name"
    message='''
python Infercnv_convert_anno_input_down.py --input 0B9YER_10x_counts.anno.raw.txt --downsample all > 0B9YER_10x_counts.anno.all.txt

--input: input anno file
--downsample: numbers of sample to include, all will include all sample, a number (50) will include a subset of samples

    '''
    print message


#050Nuc_nuclei_C01_R29   050Nuc
#050Nuc_nuclei_C01_R30   050Nuc
#050Nuc_nuclei_C01_R32   050Nuc
def parse_file(infile, downsample, samples):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
       for line in filehd:
           line = line.rstrip()
           unit = re.split(r'\t', line)
           if not line == 'x':
               #array   = re.split(r'_', unit[0])
               #sample  = array[0]
               #anno    = '%s_M' %(sample)
               
               anno = unit[1]
               if unit[1].startswith('Immune') or unit[1].startswith('Fibroblast') or unit[1].startswith('HMEC') or unit[1].startswith('GTE'):
                   #if downsample == 'all':
                   #    print '%s\t%s' %(unit[0], unit[1])
                   #elif data[unit[1]] < int(downsample):
                   #use few than 500 control
                   if data[unit[1]] < 500:
                       data[unit[1]] += 1
                       print '%s\t%s' %(unit[0], unit[1])
               elif unit[1] in ["Neg Ctrl", "Pos Ctrl", "NA"]:
                       continue
               else:
                   #if sample is specified and skip these who are not in the list
                   if len(samples.keys()) > 0:
                       if not samples.has_key(anno):
                           continue
                   if downsample == 'all':
                       print '%s\t%s' %(unit[0], anno)
                   elif data[anno] < int(downsample):
                   #sample 50 for each cluster for infercnv, too many cells will be super slow 
                       data[anno] += 1
                       print '%s\t%s' %(unit[0], anno)
    return data


#Sample  Timepoints
#1529755 Screening
#1530676 C1D14
#1541581 EOT
def read_samples(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Sample'):
                unit = re.split(r'\t',line)
                data[unit[1]] = unit[0]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-d', '--downsample')
    parser.add_argument('-s', '--sample')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0 and len(args.downsample) > 0
    except:
        usage()
        sys.exit(2)

    samples = defaultdict(lambda : str())
    if args.sample:
        samples = read_samples(args.sample)

    parse_file(args.input, args.downsample, samples)

if __name__ == '__main__':
    main()

