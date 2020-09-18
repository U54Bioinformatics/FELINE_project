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
python scRNA_Seurat_05individual_from_individual_merge_matrix.py --list FEL011027_responder_vs_non_responder.list

--list: list of patients that need to be merged

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

def read_file_pd(infile):
    file_df = pd.read_csv(infile, sep="\t", header=0) 
    return file_df

def convert_patient_id(samples):
    patients = []
    days     = []
    for sample in samples:
        patient, timepoint = re.split(r'_', sample)
        day = ''
        if timepoint == 'S':
            day = "D0"
        elif timepoint == 'M':
            day = "D14"
        elif timepoint == 'E':
            day = "D180"
        elif timepoint == 'G':
            day = "Germline"
        patients.append(patient)
        days.append(day)
    return patients, days

def get_file_df(patients, filetype):
    df_list   = []
    for p in patients:
        batch = ''
        if p == "FELINE_WES_batch1":
            batch = 'Utah'
        if p == "FELINE_WES_batch2":
            batch = 'Fulgent'
        elif p == "FELINE_WES_batch3":
            batch = 'Fulgent'
        filename = '%s%s' %(p, filetype)
        df_file = read_file_pd(filename)
        #patients, days, batches = convert_patient_id(list(df_file['Sample']), batch)
        #df_file.insert(loc=1,column="Patient", value=patients)
        #df_file.insert(loc=2,column="Timepoint", value=days)
        df_list.append(df_file)
    return df_list

#FEL027_M_TATATCCTCAAGGTGG
def read_patient_list(infile):
    file_df  = pd.read_table(infile, header=None, names=['Patient'])
    patients = list(file_df['Patient'])
    print patients
    return patients

def txt2xlsx(infile):
    file_df = pd.read_csv(infile, sep="\t", header=0)
    outxlsx = '%s.xlsx' %(infile)
    s = re.compile(r'.txt$')
    if s.search(infile):
        outxlsx = re.sub(r'.txt', r'.xlsx', infile)
    file_df.to_excel(outxlsx, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--list')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.list) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = "test_merge"

    patients = read_patient_list(args.list)
    #header      = ["Sample","TOTAL_READS","PCT_PF_READS","PCT_PF_READS_ALIGNED","MEAN_READ_LENGTH","Mean Coverage","MEAN_TARGET_COVERAGE","PCT_OFF_BAIT","PCT_USABLE_BASES_ON_TARGET","MEAN_INSERT_SIZE","STANDARD_DEVIATION","PAIR_ORIENTATION","% Duplicated"]
    #header_keep = ["Sample","TOTAL_READS","PCT_PF_READS_ALIGNED","MEAN_READ_LENGTH","MEAN_TARGET_COVERAGE","MEAN_INSERT_SIZE","STANDARD_DEVIATION","% Duplicated"]
    filetypes = ['.txt']
    for ft in filetypes:
        outfile = '%s.%s' %(args.list, ft)
        if re.compile(r'.list$').search(args.list):
            outfile = re.sub(r'.list', r'%s' %(ft), args.list)
        df_list = get_file_df(patients, ft)
        df_merged = mergefiles(df_list, on="Sample")
        patients, days = convert_patient_id(list(df_merged['Sample']))
        df_merged.insert(loc=1,column="Patient", value=patients)
        df_merged.insert(loc=2,column="Timepoint", value=days)
        df_merged = df_merged.sort_values(by=["Patient", "Timepoint"])
        df_merged.to_csv(outfile, sep="\t", index=False)
        txt2xlsx(outfile)

if __name__ == '__main__':
    main()

