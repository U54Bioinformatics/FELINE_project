#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import re
import os
import argparse

def usage():
    test="name"
    message='''
python Somatic_mutation_patient_maf.py --maf FELINE.vep.maf
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

##version 2.4
#Hugo_Symbol     Entrez_Gene_Id  Center  NCBI_Build      Chromosome      Start_Position  End_Position    Strand  Var
#FAM41C  0       .       GRCh37  1       809319  809319  +       5'Flank SNP     C       C       T       novel   
#C1orf222        0       .       GRCh37  1       1904429 1904429 +       RNA     DEL     G       G       -       nov
#PER3    0       .       GRCh37  1       7890024 7890024 +       Missense_Mutation       SNP     T       T       G
def read_maf(infile):
    data = defaultdict(lambda : str())
    maf_all        = '%s.patient.maf' %(re.sub(r'.maf', r'', infile))
    ofile_all        = open(maf_all, 'w')
    linen= 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            linen += 1
            unit = re.split(r'\t', line)
            if linen > 2:
                unit = re.split(r'\t',line)
                snp  = '%s_%s' %(unit[4], unit[5])
                unit[15] = re.split(r'_', unit[15])[0]
                newline  = '\t'.join(unit)
                print >> ofile_all, newline
            else:
                print >> ofile_all, line
    ofile_all.close()

#Chromosome      Start_Position  Reference_Allele        Tumor_Seq_Allele1       Tumor_Seq_Allele2	Tumor_Sampl
#e_Barcode	Matched_Norm_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2	t_depth	t_ref_count
#	t_alt_count	n_depth	n_ref_count	n_alt_count
#1       809319  C       C       T       FEL013_S        FEL013_G        C       -       64      57      7       22
#1       809319  C       C       T       FEL013_E        FEL013_G        C       -       68      60      8       29
def read_vep_tsv(infile):
    mutation = defaultdict(lambda : str())
    tsv_acquired   = '%s.acquired.tsv' %(re.sub(r'.tsv', r'', infile))
    tsv_maintained = '%s.maintained.tsv' %(re.sub(r'.tsv', r'', infile))
    tsv_lost       = '%s.lost.tsv' %(re.sub(r'.tsv', r'', infile))
    data_last = defaultdict(lambda : list()) 
    file_df = pd.read_table(infile, header=0)
    ofile_acquired   = open(tsv_acquired, 'w')
    ofile_maintained = open(tsv_maintained, 'w')
    ofile_lost       = open(tsv_lost, 'w')
    for i in range(file_df.shape[0]):
        snp  = '%s_%s' %(file_df['Chromosome'][i], file_df['Start_Position'][i])
        patient, sample = re.split(r'_', file_df['Tumor_Sample_Barcode'][i])
        line = '\t'.join(map(str, file_df.iloc[i,:]))
        #print line
        vaf  = 'NA'
        if not float(file_df['t_depth'][i]) == 0.0:
            vaf  = float(file_df['t_alt_count'][i])/float(file_df['t_depth'][i])
        if sample == 'S':
            data_last[snp] = [patient, sample, vaf, line]
        elif sample == 'E':
            if data_last[snp][2] == 'NA' or vaf == 'NA':
                continue
            if float(data_last[snp][2]) <= 0.10 and float(vaf) >= 0.20:
                # acquired
                #print >> ofile_acquired, data_last[snp][3]
                print >> ofile_acquired, line
                mutation[snp] = 'acquired'
            elif float(data_last[snp][2]) >= 0.10 and float(vaf) >= 0.10:
                # maintained
                print >> ofile_maintained, data_last[snp][3]
                #print >> ofile_maintained, line
                mutation[snp] = 'maintained'
            elif float(data_last[snp][2]) >= 0.20 and float(vaf) <= 0.10:
                # lost
                print >> ofile_lost, data_last[snp][3]
                #print >> ofile_lost, line
                mutation[snp] = 'lost'
    ofile_acquired.close()
    ofile_maintained.close()
    ofile_lost.close()
    return mutation

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
    parser.add_argument('--maf')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.maf) > 0
    except:
        usage()
        sys.exit(2)

    read_maf(args.maf)

if __name__ == '__main__':
    main()

