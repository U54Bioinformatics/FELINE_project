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

#FELINE_patient_1_46.clinical_data.txt
#Sample  orig.ident      Patient.Study.ID        ARM     Treatment       Ribo    Ki67_Response   Response_v2     R
#FEL001P120      FEL001  001-120 B       letrozole + ribo        1       Sensitive       Non-responder   Progressive
#FEL002P106      FEL002  001-106 C       letrozole + ribo        1       Sensitive       Responder       Partial re
def parse_clinic_data(infile):
    data = defaultdict(lambda : list())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        data[file_df['orig.ident'][i]] = [file_df['ARM'][i], file_df['Response'][i]]
    return data

#Sample	Purity	Ploidy	Ploidy_Int
#FEL013_E	0.493025392515	3.14281122457	4
#FEL013_S	0.746287834318	3.56908668798	4
def parse_ploidy(infile):
    data = defaultdict(lambda : list())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        data[file_df['Sample'][i]] = [file_df['Purity'][i], file_df['Ploidy_Int'][i]]
    return data

#FACETS output format:
#data frame of segment summaries pre and post clustering of segments. The columns are: chrom the chromosome to which the segment belongs; seg the segment number; num.mark the number of SNPs in the segment; nhet the number of SNPs that are deemed heterozygous; cnlr.median the median log-ratio of the segment; mafR the log-odds-ratio summary for the segment; segclust the segment cluster to which segment belongs; cnlr.median.clust the median log-ratio of the segment cluster; mafR.clust the log-odds-ratio summary for the segment cluster; cf the cellular fraction of the segment; tcn the total copy number of the segment; lcn the minor copy number of the segment.
#chrom  start   end     seg     num.mark        nhet    cnlr.median     mafR    segclust        cnlr.median.clust       mafR.clust      cf.em   tcn.em  lcn.em
#output:
#chromosome     start   end     copy_number     major_cn        minor_cn        cellular_prevalence
#1      1       100000  2       1       1       0.81
#1      141500000       148899999       2       1       1       0.81
#1      148900000       163507499       2       1       1       0.81
def parse_facets_cnv(infile, outfile, ploidy):
    data = defaultdict(lambda : str())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    #print file_df["cf.em"][2]
    file_df = file_df[["chrom", "start", "end", "nhet", "cnlr.median", "tcn.em", "lcn.em", "cf.em"]]
    file_df.columns = ["chromosome", "start", "end", "nhet", "cnlr.median", "copy_number", "minor_cn", "cellular_prevalence"]
    copy_call = []
    minor_cn  = []
    cp        = []
    for i in range(file_df.shape[0]):
        #print file_df['start'][i]
        try:
            if pd.isnull(file_df['minor_cn'][i]):
                #print "flag1 %s" %(file_df['minor_cn'][i])
                minor_cn.append("NA")
            else:
                minor_cn.append(file_df['minor_cn'][i])
        except:
            #print "flag3 %s" %(file_df['minor_cn'][i])
            minor_cn.append("NA")
        try:
            if pd.isnull(file_df['cellular_prevalence'][i]):
                cp.append("NA")
            else:
                cp.append(file_df['cellular_prevalence'][i])
        except:
            cp.append("NA")
        
        if float(file_df['cnlr.median'][i]) >= 0.2 and int(file_df['copy_number'][i]) > int(ploidy[1]):
            # total copy number gain
            copy_call.append("Gain")
        elif float(file_df['cnlr.median'][i]) <= -0.2 and int(file_df['copy_number'][i]) < int(ploidy[1]):
            # total copy number loss
            if float(file_df['minor_cn'][i]) == 0.0:
                copy_call.append("LOH")
            else:
                copy_call.append("Loss")
        else:
            # total copy number normal
            if file_df['minor_cn'][i] != "NA":
                #if float(file_df['minor_cn'][i]) == 0.0 and int(ploidy[1]) == int(file_df['copy_number'][i]):
                if float(file_df['minor_cn'][i]) == 0.0 and int(file_df['copy_number'][i]) > 0:
                    copy_call.append("LOH")
                else: 
                    copy_call.append("Normal")
            else:
                copy_call.append("Normal")
    file_df['CNV_call'] = copy_call
    file_df['minor_cn'] = minor_cn 
    file_df['cellular_prevalence'] = cp
    #print file_df
    file_df = file_df[["chromosome", "start", "end", "CNV_call", "nhet", "cnlr.median", "copy_number", "minor_cn", "cellular_prevalence"]]
    #file_df = file_df.dropna()
    file_df = file_df.astype({"start": "int64"})
    file_df = file_df.astype({"end": "int64"})
    file_df.to_csv(outfile, sep="\t", index=False, header=False)


def prepare_cnv_driver(infile, driver_file1, driver_file2, drivers):
    data = defaultdict(lambda : str())
    file_df = pd.read_csv(infile, header=None, sep="\t")
    file_df.columns = ["chromosome", "start", "end", "CNV_call", "nhet", "cnlr.median", "copy_number", "minor_cn", "cellular_prevalence", "Chrom2", "Start2", "End2", "Gene.ID"]
    ofile1   = open(driver_file1, 'w')
    ofile2   = open(driver_file2, 'w')
    header   = "gene\tchromosome\tstart\tend\tCNV_call\tnhet\tcnlr.median\ttotal_cn\tminor_cn\tcellular_prevalence"
    print >> ofile1, header
    print >> ofile2, header
    outlines = defaultdict(lambda : str())
    for i in range(file_df.shape[0]):
        gene  = file_df["Gene.ID"][i]
        #if drivers.has_key(gene) and file_df["CNV_call"][i] != "Normal":
        if drivers.has_key(gene):
            out_line = [gene, file_df["chromosome"][i], file_df["start"][i], file_df["end"][i], file_df["CNV_call"][i], file_df["nhet"][i], file_df["cnlr.median"][i], file_df["copy_number"][i], file_df["minor_cn"][i], file_df["cellular_prevalence"][i]]
            print >> ofile1, '\t'.join(map(str, out_line))
            if gene in ["AKT1", "AKT3", "ERBB4", "PIK3CA", "ESR1", "CDK6", "FGFR1", "FGFR2", "CCNE2", "CCND1", "MYC", "RB1", "TP53", "ERBB2", "NF1"]:
                #print >> ofile2, '\t'.join(map(str, out_line))
                outlines[gene] = '\t'.join(map(str, out_line))
    ofile1.close()
    for g in ["AKT1", "AKT3", "ERBB4", "PIK3CA", "ESR1", "CDK6", "FGFR1", "FGFR2", "CCNE2", "CCND1", "MYC", "RB1", "TP53", "ERBB2", "NF1"]:
        print >> ofile2, outlines[g]
    ofile2.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--facets')
    parser.add_argument('--sample')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.facets) > 0
    except:
        usage()
        sys.exit(2)

    ploidy = parse_ploidy("FELINE_FACETS_tumor_purity.txt")

    clinic = parse_clinic_data("FELINE_patient_1_46.clinical_data.txt")

    # bed file facets results
    bed_cnv = '%s.facets.seg_cnv_call.bed' %(args.sample)
    if args.facets:
        parse_facets_cnv(args.facets, bed_cnv, ploidy[args.sample])


    # hg19 gene bed
    bed_hg19 = 'hg19.RefSeq.NM_pos_unique_sort.bed'
    # bedtoools overlapping
    overlap_file = '%s.facets.seg_cnv_call.overlap' %(args.sample)
    cmd = 'bedtools window -w 1 -a %s -b %s > %s' %(bed_cnv, bed_hg19, overlap_file)
    os.system(cmd)
    # prepare driver file
    gene_list = "Cancer_genes.drivers_resistent.list"
    drivers   = read_driver_file(gene_list)

    # write driver gene cnv output
    driver_file1 = '%s.facets.gene_cnv_call.txt' %(args.sample)
    driver_file2 = '%s.facets.gene_cnv_call.short_list.txt' %(args.sample)
    prepare_cnv_driver(overlap_file, driver_file1, driver_file2, drivers)
  
if __name__ == '__main__':
    main()

