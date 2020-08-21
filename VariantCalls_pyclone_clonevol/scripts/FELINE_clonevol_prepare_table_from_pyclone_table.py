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
python CircosConf.py --input circos.config --output pipe.conf
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

#Gene Symbol
#ABL1
#ACVR1B
#AKT1
def read_driver_file(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        data[file_df['Gene Symbol'][i]] = '1'
    return data

#mutation_id     sample_id       cluster_id      cellular_prevalence     cellular_prevalence_std variant_allel
#10_100177815    FEL013_S        0       0.20871286667966688     0.05276455573816322     0.10416666666666667
#10_100177815    FEL013_E        0       0.032531644137520396    0.08092668011152085     0.0
def read_frequency_file_ccf(infile):
    data = defaultdict(lambda : defaultdict(lambda : str()))
    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        # mutation -> sample -> frequency
        #print file_df['mutation'][i], file_df['FEL034_E'][i], file_df['FEL034_S'][i]
        data[file_df['mutation_id'][i]][file_df['sample_id'][i]] = '%.2f' %(100 * float(file_df['cellular_prevalence'][i]))
        data[file_df['mutation_id'][i]]["cluster_id"] = file_df['cluster_id'][i]
    return data



#cancer_cell_frequency.txt
#mutation        FEL034_E        FEL034_S
#2_91894137      0.440490847951  0.00445159552673
#2_136742966     0.0103726051662 0.30874531754
#variant_allele_frequency.txt
#mutation        FEL034_E        FEL034_S
#10_132891385    0.140350877193  0
#10_29598895     0.0105263157895 0.116071428571
#10_74631403     0.147540983607  0.216216216216
def read_frequency_file(infile):
    data = defaultdict(lambda : defaultdict(lambda : str()))
    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        # mutation -> sample -> frequency
        #print file_df['mutation'][i], file_df['FEL034_E'][i], file_df['FEL034_S'][i]
        for j in range(1, len(file_df.columns)):
            #print file_df['mutation'][i], file_df.columns[j], file_df[file_df.columns[j]][i]
            data[file_df['mutation'][i]][file_df.columns[j]] = '%.2f' %(100 * float(file_df[file_df.columns[j]][i]))
    return data, file_df.columns[1:]



#cluster.frequency.final.txt
#ClusterID	Cluster	Num Mutations	FEL034_E	FEL034_S
#1	0	40	0.0103726061061	0.308745354414
#2	1	36	0.440490782261	0.00445159571245
#3	2	29	0.943601191044	0.966029644012
#4	3	53	0.287135064602	0.246382638812
def read_cluster_file(infile):
    data = defaultdict(lambda : str())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        #print file_df[0][i]
        # original pyclone cluster ID -> clonevol cluster ID
        data[file_df['Cluster'][i]] = file_df['ClusterID'][i]
    return data

#mutations.txt
#Mutation        Cluster Gene ID Gene Symbol
#10_29598895     0       ENSG00000120563,ENST00000375500,ENST00000494304 LYZL1
#10_74631403     3       ENST00000605416 MCU
#10_88822682     0       ENST00000474574 GLUD1
def process_mutation_file(infile, sample, cluster, vaf, ccf, driver, prefix):
    # ccf outfile
    outfile_ccf = '%s.clonevol_ccf.txt' %(prefix)
    ofile_ccf = open(outfile_ccf, 'w')
    sample_ccf= '\t'.join(sample)
    sample_vaf= '\t'.join(['%s.vaf' %(x) for x in sample])
    print >> ofile_ccf, 'mutation\tcluster\tgene\tis.driver\t%s\t%s' %(sample_ccf, sample_vaf)  
    # vaf outfile
    outfile_vaf = '%s.clonevol_vaf.txt' %(prefix)
    ofile_vaf = open(outfile_vaf, 'w')
    sample_ccf= '\t'.join(['%s.ccf' %(x) for x in sample])
    sample_vaf= '\t'.join(sample)
    print >> ofile_vaf, 'mutation\tcluster\tgene\tis.driver\t%s\t%s' %(sample_ccf, sample_vaf)   

    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        #print file_df[0][i]
        # original pyclone cluster ID -> clonevol cluster ID
        #data[file_df['Cluster'][i]] = file_df['ClusterID'][i]
        if ccf.has_key(file_df['Mutation'][i]):
            # pyclone id
            cluster_id_raw = ccf[file_df['Mutation'][i]]["cluster_id"] 
            # assigned cluster for clonevol, starting from 1 and continuous
            # skip these clusters that are not in final cluster tables
            if not cluster.has_key(cluster_id_raw):
                continue
            cluster_id = cluster[cluster_id_raw]
            gene       = file_df['Gene Symbol'][i]
            isdriver     = 'TRUE' if driver.has_key(gene) else 'FALSE' 
            if ccf.has_key(file_df['Mutation'][i]) and vaf.has_key(file_df['Mutation'][i]):
                ccf_list = []
                vaf_list = []
                for s in sample:
                    ccf_list.append(ccf[file_df['Mutation'][i]][s])
                    vaf_list.append(vaf[file_df['Mutation'][i]][s])
                ccf_line = '\t'.join(ccf_list)
                vaf_line = '\t'.join(vaf_list)
                newline = '%s\t%s\t%s\t%s\t%s\t%s' %(file_df['Mutation'][i], cluster_id, gene, isdriver, ccf_line, vaf_line)
                print >> ofile_ccf, newline
                print >> ofile_vaf, newline
    ofile_ccf.close()
    ofile_vaf.close()

def read_large_matrix(infile):
    # determine and optimize dtype
    # Sample 100 rows of data to determine dtypes.
    file_test = pd.read_csv(infile, sep="\t", header=0, nrows=100)
    float_cols = [c for c in file_test if file_test[c].dtype == "float64"]
    int_cols = [c for c in file_test if file_test[c].dtype == "int64"]
    if float_cols > 0:
        dtype_cols = {c: np.float32 for c in float_cols}
    elif int_cols > 0:
        dtype_cols = {c: np.int32 for c in int_cols}
    file_df = pd.read_csv(infile, sep="\t", header=0, dtype=dtype_cols)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pyclone_dir')
    parser.add_argument('--driver')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.pyclone_dir) > 0
    except:
        usage()
        sys.exit(2)

    # name output prefix
    if not args.output:
        temp  = args.pyclone_dir
        args.output = temp.rstrip('/')

    driver_file   = args.driver
    cluster_file  = '%s/tables/pyclone.table_cluster.final.txt' %(args.pyclone_dir)
    mutation_file = '%s/mutations.txt' %(args.pyclone_dir)
    ccf_file      = '%s/tables/pyclone.table_loci.txt' %(args.pyclone_dir)
    vaf_file      = '%s/data_files/variant_allele_frequency.txt' %(args.pyclone_dir)

    driver      = read_driver_file(driver_file)
    ccf         = read_frequency_file_ccf(ccf_file)
    vaf, sample = read_frequency_file(vaf_file)
    cluster     = read_cluster_file(cluster_file)
    process_mutation_file(mutation_file, sample, cluster, vaf, ccf, driver, args.output)

if __name__ == '__main__':
    main()

