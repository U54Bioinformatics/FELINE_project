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
python FELINE_Step5_merge_fishplot.py --input FELINE_patient_info.sorted.txt

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

#Pateint cnv
#FEL013  FACETS
#FEL014  sequenza
#FEL015  sequenza
def read_pyclone_model(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        data[file_df['Pateint'][i]] = file_df['cnv'][i]
    return data

#PID     Sample  orig.ident      ARM     Response        Cell    Percet  Day0    Day14   Day180  Pre     Post
#P15     FEL015P105      FEL015  A       Non-responder   972     44.73   43      83      846     Yes     Yes
#P20     FEL020P115      FEL020  A       Non-responder   9677    97.87   3893    3036    2748    Yes     Yes
def read_clinic(infile, prefix, response, model):
    data = defaultdict(lambda : str())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    trees   = []
    tree_cmd= []
    labels  = []
    for i in range(file_df.shape[0]):
        sample = file_df['orig.ident'][i]
        arm    = file_df['ARM'][i]
        #reponse= file_df['Response'][i]
        response = "NR" if file_df['Response'][i] == "Non-responder" else "R"
        if model.has_key(sample) and not sample in ["FEL021", "FEL029"]:
            cnv  = model[sample]
            cmd1 = 'p="%s"' %(sample)
            cmd2 = 'cnv="%s"' %(cnv)
            cmd3 = 'tree%s  = ggdraw() + draw_image(paste0(p, "_pyclone_analysis_",cnv ,"_20_5_0.05_1_10000_none.clonevol_ccf.clonevol_figure4_fishplot.png"))'  %(i)
            tree_cmd.append(cmd1)
            tree_cmd.append(cmd2)
            tree_cmd.append(cmd3)
            trees.append("tree%s" %(i))
            lab = re.sub(r'FEL0', r'P', sample)
            labels.append('"%s (Arm %s, %s)"' %(lab, arm, response))
            #labels.append('"%s (Arm %s)"' %(lab, arm))

    head='''library(ggplot2)
library(cowplot)
library(magick)
'''
    end='''
pdf("%s.pdf", height=22, width=16)
plot_grid(%s, nrow=6, ncol=4, labels=c(%s), label_size=20)
dev.off()
''' %(prefix, ",".join(trees), ','.join(labels))

    ofile = open("%s.R" %(prefix), 'w')
    print >> ofile, head
    for cmd in tree_cmd:
        print >> ofile, cmd
    print >> ofile, end
    ofile.close()  

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--pyclone_cnv')
    parser.add_argument('--response')
    parser.add_argument('--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = "FELINE_fishplot"
   
    if not args.response:
        args.response = "Non-resoponder"
 
    if not args.pyclone_cnv:
        args.pyclone_cnv = "FELINE_pyclone_model.txt"    

    # read pyclone model
    model = read_pyclone_model(args.pyclone_cnv)

    # write merge script
    read_clinic(args.input, args.project, args.response, model)

if __name__ == '__main__':
    main()

