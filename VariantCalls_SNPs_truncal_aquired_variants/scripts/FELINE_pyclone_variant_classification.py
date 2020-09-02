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
python FELINE_pyclone_variant_classification.py --patient FEL013 --driver Cancer_genes.ID.list --cnv sequenza

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

#get gene annotation for variants
#FELINE_variants/FEL013_variant_calls_04_purity40_totalRD20_minread5_minVAF0.05.common_snp_info.txt
#Chrom   Pos     Ref     Alt     FEL013_E        FEL013_S        Func.refGene    Gene.refGene
#1       985450  G       A       165/1/0.006     131/7/0.050     intronic        AGRN                            
#1       1904429 G       A       193/0/0.000     148/51/0.256    exonic  CFAP74          synonymous SNV        
def get_variant_gene(infile):
    #infile = "FELINE_variants/%s_variant_calls_04_purity40_totalRD20_minread5_minVAF0.05.common_snp_info.txt" %(patient)
    data = defaultdict(lambda : str())
    large_impact = defaultdict(lambda : list())
    file_df = pd.read_csv(infile, header=0, sep="\t")
    #print file_df.shape
    for i in range(file_df.shape[0]):
        gene = file_df["Gene.refGene"][i]
        variant = '%s_%s' %(file_df["Chrom"][i], file_df["Pos"][i])
        data[variant] = gene
        #if file_df["Func.refGene"][i] in ["exonic;splicing", "exonic", "splicing", "ncRNA_exonic", "ncRNA_splicing"]:
        if file_df["Func.refGene"][i] in ["exonic;splicing", "exonic", "splicing"]:
            if not file_df["ExonicFunc.refGene"][i] in ["synonymous SNV", "unknown"]:
                large_impact[variant] = [file_df["Func.refGene"][i], file_df["ExonicFunc.refGene"][i]]
    return data, large_impact


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
        data[str(file_df['ClusterID'][i])] = [int(100*file_df[file_df.columns[3]][i]), int(100*file_df[file_df.columns[4]][i])]
        #print data[str(file_df['ClusterID'][i])]
    return data, file_df.columns[3:]

#mutation        cluster gene    is.driver       FEL013_E        FEL013_S        FEL013_E.vaf    FEL013_S.vaf
#10_100177815    2       HPS1    FALSE   3.25    20.87   0.00    10.42
#10_105037095    2       INA     FALSE   0.73    20.64   0.00    7.55
#10_112838100    4       ADRA2A  FALSE   14.36   17.82   2.31    7.35
#10_112838101    2       ADRA2A  FALSE   0.84    19.96   0.52    6.47
#samples: [FEL045, FEL045]
#impact: large_impact[variant] = [file_df["Func.refGene"][i], file_df["ExonicFunc.refGene"][i]]
#cluster: cluster_id -> [cluster ccf of FEL045_M, cluster ccf of FEL045_S]
def process_ccf_file(infile, patient, cluster, gene, impact, samples, driver, outfile):
    ofile = open(outfile, 'w')
    ccf_temp = ['%s.ccf' %(x) for x in samples]
    vaf_temp = ['%s.vaf' %(x) for x in samples]
    ccf_header = ['%s' %(x) for x in samples]
    vaf_header = ['%s.vaf' %(x) for x in samples]
    print >> ofile, 'patient\tmutation\tgene\ttype\tcluster\tcluster_ccf_post\tcluster_ccf_pre\t%s\t%s\tdriver\tFunc.refGene\tExonicFunc.refGene' %('\t'.join(ccf_temp), '\t'.join(vaf_temp))
    file_df = pd.read_csv(infile, header=0, sep="\t")
    for i in range(file_df.shape[0]):
        mutation_type = 'NA'
        cluster_id    = str(file_df['cluster'][i])
        if cluster.has_key(cluster_id):
            if int(file_df['cluster'][i]) == 1:
                # truncal mutations, need to have ccf>=70% at both timepoints
                if int(cluster[cluster_id][0]) >= 70 and int(cluster[cluster_id][1]) >= 70:
                    mutation_type = "truncal"
                if float(file_df[ccf_header[1]][i]) >= 70 and float(file_df[ccf_header[0]][i]) >= 70:
                    mutation_type = "truncal"
            else:
                # acquired mutations, need to have ccf<=10% at day 0 and >= 20% at day 14/180 
                if int(cluster[cluster_id][1]) <= 10 and int(cluster[cluster_id][0]) >= 20:
                    mutation_type = "acquired"
                if float(file_df[ccf_header[1]][i]) <= 10 and float(file_df[ccf_header[0]][i]) >= 20:
                    mutation_type = "acquired" 
        #print cluster_id, mutation_type
        large_impact = 0
        if impact.has_key(file_df['mutation'][i]):
            large_impact = 1
        if large_impact == 1:
            # only output large impact mutations
            #ccf_header = ['%s' %(x) for x in samples]
            #vaf_header = ['%s.vaf' %(x) for x in samples]
            driver_flag= "FALSE"
            if driver.has_key(gene[file_df['mutation'][i]]):
                driver_flag = "TRUE"
            output_list= [patient, file_df['mutation'][i], gene[file_df['mutation'][i]], mutation_type, cluster_id, cluster[cluster_id][0], cluster[cluster_id][1],  file_df[ccf_header[0]][i], file_df[ccf_header[1]][i], file_df[vaf_header[0]][i], file_df[vaf_header[1]][i], driver_flag, impact[file_df['mutation'][i]][0], impact[file_df['mutation'][i]][1]]
            print >> ofile, '\t'.join(map(str, output_list))
    ofile.close()

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
    parser.add_argument('--patient')
    parser.add_argument('--pyclone_dir')
    parser.add_argument('--mutation_dir')
    parser.add_argument('--cnv')
    parser.add_argument('--driver')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.patient) > 0
    except:
        usage()
        sys.exit(2)

    # name output prefix
    #if not args.output:
    #    temp  = args.pyclone_dir
    #    args.output = temp.rstrip('/')

    if not args.pyclone_dir:
        args.pyclone_dir = "./VariantCalls_pyclone_clonevol"

    if not args.mutation_dir:
        args.mutation_dir = "./FELINE_variants"

    driver_file   = args.driver
    patient_folder= "%s_pyclone_analysis_%s_20_5_0.05_1_10000_none" %(args.patient, args.cnv)
    cluster_file  = '%s/%s/tables/pyclone.table_cluster.final.txt' %(args.pyclone_dir, patient_folder)
    mutation_file = '%s/%s_variant_calls_04_purity40_totalRD20_minread5_minVAF0.05.common_snp_info.txt' %(args.mutation_dir, args.patient)
    ccf_file      = '%s/%s.clonevol_ccf.txt' %(args.pyclone_dir, patient_folder)

    driver          = read_driver_file(driver_file)
    cluster, sample = read_cluster_file(cluster_file)
    gene, impact    = get_variant_gene(mutation_file)

    #outfile = '%s.variant_info.%s.txt' %(args.patient, args.cnv)
    outfile = '%s.variant_info.txt' %(args.patient)
    process_ccf_file(ccf_file, args.patient, cluster, gene, impact, sample, driver, outfile)

if __name__ == '__main__':
    main()

