#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import pandas
import re
import os
import argparse
import glob

def usage():
    test="name"
    message='''
python Run_merge_tree_helper.py --input FELINE_patient_info.sorted.txt

    '''
    print message


def mergefiles(dfs=[], on=''):
    """Merge a list of files based on one column"""
    if len(dfs) == 1:
         return "List only have one element."

    elif len(dfs) == 2:
        df1 = dfs[0]
        df2 = dfs[1]
        df = df1.merge(df2, on=on)
        return df

    # Merge the first and second datafranes into new dataframe
    df1 = dfs[0]
    df2 = dfs[1]
    df = dfs[0].merge(dfs[1], on=on)

    # Create new list with merged dataframe
    dfl = []
    dfl.append(df)

    # Join lists
    dfl = dfl + dfs[2:] 
    dfm = mergefiles(dfl, on)
    return dfm

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

#PID     Sample  orig.ident      ARM     Response        Cell    Percet  Day0    Day14   Day180  Pre     Post
#P15     FEL015P105      FEL015  A       Non-responder   972     44.73   43      83      846     Yes     Yes
#P20     FEL020P115      FEL020  A       Non-responder   9677    97.87   3893    3036    2748    Yes     Yes
def read_clinic(infile, prefix, response):
    data = defaultdict(lambda : str())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    trees   = []
    tree_cmd= []
    for i in range(file_df.shape[0]):
        sample = file_df['Sample'][i]
        if file_df['Response'][i] == response and sample != "FEL032P145":
            cmd1 = 'p="%s"' %(sample)
            cmd2 = 'tree%s  = ggdraw() + draw_image(paste0("./MEDALT_tree_png/", p, ".tree.png"))'  %(i)
            tree_cmd.append(cmd1)
            tree_cmd.append(cmd2)
            trees.append("tree%s" %(i))

    head='''library(ggplot2)
library(cowplot)
library(magick)
'''
    if response == "Non-responder":
        end='''
pdf("%s.%s.pdf", height=22, width=10)
plot_grid(%s, nrow=6, ncol=2)
dev.off()
''' %(prefix, response, ",".join(trees))
    elif response == "Responder":
        end='''
pdf("%s.%s.pdf", height=22, width=10)
plot_grid(%s, nrow=6, ncol=2)
dev.off()
''' %(prefix, response, ",".join(trees))
  

    ofile = open("%s.%s.R" %(prefix, response), 'w')
    print >> ofile, head
    for cmd in tree_cmd:
        print >> ofile, cmd
    print >> ofile, end
    ofile.close()  

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-r', '--response')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = "FELINE_scRNA_tree_25patients"
   
    if not args.response:
        args.response = "Non-resoponder"

    read_clinic(args.input, args.project, args.response)

if __name__ == '__main__':
    main()

