#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import re
import os
import argparse
import glob

def usage():
    test="name"
    message='''
python Merge_multiple_matrix.py --input ./test/ --output test_merge

--input: folder with individual file to be merged
--output: prefix of merged file

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

#FELINE_patient_1_46.clinical_data.txt
#Sample  orig.ident      Patient.Study.ID        ARM     Treatment       Ribo    Ki67_Response   Response_v2     R
#FEL001P120      FEL001  001-120 B       letrozole + ribo        1       Sensitive       Non-responder   Progressive
#FEL002P106      FEL002  001-106 C       letrozole + ribo        1       Sensitive       Responder       Partial re
def parse_clinic_data(infile):
    data = defaultdict(lambda : list())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        for tp in ["S", "M", "E"]:
            data["%s_%s" %(file_df['orig.ident'][i], tp)] = [file_df['ARM'][i], file_df['Response'][i]]
    return data

#Sample Purity  Ploidy  Ploidy_Int
#FEL013_E       0.493025392515  3.14281122457   4
#FEL013_S       0.746287834318  3.56908668798   4
def parse_ploidy(infile):
    data = defaultdict(lambda : list())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        data[file_df['Sample'][i]] = [file_df['Purity'][i], file_df['Ploidy_Int'][i]]
    return data

def convert_sample(sample):
    sample_new = re.sub(r'FEL0', r'P', sample)
    if "_S" in sample_new:
        sample_new = re.sub(r'_S', r'_T1', sample_new)
    if "_M" in sample_new:
        sample_new = re.sub(r'_M', r'_T2', sample_new)
    if "_E" in sample_new:
        sample_new = re.sub(r'_E', r'_T3', sample_new)
    return sample_new
    

#FELINE_FACETS_gene_cnv/FEL014_S.facets.gene_cnv_call.short_list.txt
#gene    chromosome      start   end     CNV_call        nhet    cnlr.median     total_cn        minor_cn    
#AKT1    14      95557012        105996000       Normal  149     0.104034517043412       4       2.0     0.4406923
#AKT3    1       164529120       249212569       Normal  859     -0.0350683584215823     3       1.0     0.63523084
def read_file_pd(infile, ploidy, clinic):
    file_fn = os.path.split(infile)[1]
    sample  = re.split(r'\.', file_fn)[0]
    file_df = pd.read_csv(infile, sep="\t", header=0)
    # add a new column for patient sample
    sample_1 = re.sub(r'FEL0', r'P', sample)
    #sample_2 = convert_sample(sample)
    sample_2 = sample
    file_df["sample"]  = [sample_1 for i in range(file_df.shape[0])]
    file_df["sample1"] = [sample_2 for i in range(file_df.shape[0])]
    # add a new column for arm
    arm = clinic[sample][0]
    file_df["arm"] = [arm for i in range(file_df.shape[0])]
    # add a new column for response
    response = clinic[sample][1]
    file_df["response"] = [response for i in range(file_df.shape[0])]
    # add a new column for purity
    purity = ploidy[sample][0]
    file_df["purity"] = [purity for i in range(file_df.shape[0])]
    # add a new column for ploidy
    ploidy = ploidy[sample][1]
    file_df["ploidy"] = [ploidy for i in range(file_df.shape[0])]
    # change Normal to N
    file_df["CNV_call"] = [re.sub(r'Normal', r'N', file_df['CNV_call'][i]) for i in range(file_df.shape[0])]
    # reorder column to make sample as first column
    header  = ["sample", "sample1", "arm", "response", "purity", "ploidy"]
    header.extend(file_df.columns[:-6])
    file_df = file_df[header]
    return file_df

def txt2xlsx(infile):
    file_df = pd.read_csv(infile, sep="\t", header=0)
    outxlsx = '%s.xlsx' %(infile)
    s = re.compile(r'.txt$')
    if s.search(infile):
        outxlsx = re.sub(r'.txt', r'.xlsx', infile)
    file_df.to_excel(outxlsx, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = "test_merge"

    ploidy = parse_ploidy("FELINE_FACETS_tumor_purity.txt")

    clinic = parse_clinic_data("FELINE_patient_1_46.clinical_data.txt")


    file_list = glob.glob('%s/*.gene_cnv_call.txt' %(args.input))
    df_list   = []
    for f in sorted(file_list):
        df_file = read_file_pd(f, ploidy, clinic)
        df_list.append(df_file)
    
    # merge gene cnv from all patient sample into a big meta table (source data)
    outfile   = '%s.gene_cnv_call.txt' %(args.output)
    #df_merged = mergefiles(df_list, on="Gene Set")
    df_merged = pd.concat(df_list)
    df_merged.to_csv(outfile, sep="\t", index=False)
    txt2xlsx(outfile)

    #
    #genes   sample_name     CNV
    #AL627309.1      FEL013_E        Amp
    #CICP27  FEL013_E        Amp
    #FAM138A FEL013_E        Amp
    #MIR1302-10      FEL013_E        Amp
    outfile4maftools = '%s.gene_cnv_call4maftools.txt' %(args.output)
    df_maftools = df_merged[["gene", "sample1", "CNV_call"]] 
    df_maftools.to_csv(outfile4maftools, sep="\t", index=False)


    
    rmdup = defaultdict(lambda : int())
    df_maftools = pd.read_csv(outfile4maftools, header=0, sep="\t")
    ofile = open(outfile4maftools, 'w')
    print >> ofile, "genes\tsample_name\tCNV"
    print df_maftools.columns
    for i in range(df_maftools.shape[0]):
        if str(df_maftools['CNV_call'][i]) in ["N", "Normal"]:
            continue
        elif str(df_maftools['CNV_call'][i]) in ["Loss", "LOH"]:
            call = "%s.%s" %(df_maftools['gene'][i], df_maftools['sample1'][i])
            if not rmdup.has_key(call):
                print >> ofile, "%s\t%s\t%s" %(df_maftools['gene'][i], df_maftools['sample1'][i], "Del")
                rmdup[call] = 1
        elif str(df_maftools['CNV_call'][i]) in ["Gain"]:
            call = "%s.%s" %(df_maftools['gene'][i], df_maftools['sample1'][i])
            if not rmdup.has_key(call):
                print >> ofile, "%s\t%s\t%s" %(df_maftools['gene'][i], df_maftools['sample1'][i], "Amp")
                rmdup[call] = 1
    ofile.close()
 
if __name__ == '__main__':
    main()

