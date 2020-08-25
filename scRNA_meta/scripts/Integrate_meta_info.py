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
python Prepare_ssGSEA_data_per_pathway.py --anno FEL039P201_10x_cell_metadata.UMAPcluster.SingleR_anno_cluster.revised_anno.infercnv_predict.txt --input FEL039P201_10x.ssGSEA.scores.txt --patient FEL039P201

Split ssGESA score by pathway for each patient

Output:
Pathway Cell.ID Sample  Timepoint       Celltype        Score
ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_DN   FEL039_E_AACAAGACAGCGTATT       FEL039P201      FEL039P201_E    Epithelial cells        -0.0429953840681352
ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_DN   FEL039_E_AACGAAATCGGTTGTA       FEL039P201      FEL039P201_E    Epithelial cells        0.00915629007790477
ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_DN   FEL039_E_AAGAACAGTCGTTGGC       FEL039P201      FEL039P201_E    Cancer cells    -0.0214135210119024
ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_DN   FEL039_E_AATAGAGCATGAAGGC       FEL039P201      FEL039P201_E    Epithelial cells        -0.0107961024917405
ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_DN   FEL039_E_ACGGGTCCATCCGATA       FEL039P201      FEL039P201_E    Cancer cells    -0.0293698244736759
ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_DN   FEL039_E_AGACCCGTCAATCTCT       FEL039P201      FEL039P201_E    Cancer cells    -0.01691212921179
ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_DN   FEL039_E_AGCATCATCGTTACCC       FEL039P201      FEL039P201_E    Epithelial cells        -0.0257454404888769
ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_DN   FEL039_E_AGCGCTGGTCTCTCAC       FEL039P201      FEL039P201_E    Epithelial cells        0.00982595803571441
ACEVEDO_LIVER_CANCER_WITH_H3K27ME3_DN   FEL039_E_AGGCTGCAGGGCAGAG       FEL039P201      FEL039P201_E    Cancer cells    -0.00908561791133496

    '''
    print message


#Sample  orig.ident      Patient.Study.ID        ARM     Treatment       Ribo    Ki67_Response   Response_v2     Res
#FEL001P120      FEL001  001-120 B       letrozole + ribo        1       Sensitive       Non-responder   Progressive
#FEL002P106      FEL002  001-106 C       letrozole + ribo        1       Sensitive       Responder       Partial res

def read_clinical(infile):
    data = defaultdict(lambda : list())
    file_df = pd.read_table(infile, header=0)
    for i in range(file_df.shape[0]):
        data[file_df['orig.ident'][i]] = [file_df['ARM'][i], file_df['Treatment'][i], file_df['Ki67_Response'][i], file_df['Response_v3b'][i], file_df['Response'][i]]
    return data

#Cell.ID orig.ident      nCount_RNA      nFeature_RNA    Sample  Total.Reads     Total.Reads..Thousands. Total.Non.Mitochondria.Reads    Total.Non.Mitochondria.Reads..Thousands
#FEL034_M_GGGCCATCAGTCAGCC       FEL034  4061    2112    FEL034P601_M    4062    4.062   4055    4.055   2113    -0.00456805214465       -0.048958939064 G1      0.2
#FEL034_M_TAGGTACGTATGAGGC       FEL034  2234    1257    FEL034P601_M    2235    2.235   2221    2.221   1258    0.010606768524  0.0133580085921 G2M     0.6
def read_cell_anno(infile, clinical_data, outfile):
    data = defaultdict(lambda : str())
    file_df = pd.read_table(infile, header=0)
    header  = ["Cell.ID", 'orig.ident', 'nCount_RNA', 'nFeature_RNA', "Sample", "Timepoint", "Day", "Celltype", "Celltype_subtype", "Platform", "ARM", "Treatment", "Ki67_Response", "Response_v3b", "Response"]
    #file_df = file_df[header]
    #file_df.columns = ["Cell.ID", 'nCount_RNA', 'nFeature_RNA', "Timepoint", "Celltype", "Celltype_subtype", "Platform"]
    ofile = open(outfile, 'w')
    print >> ofile, '\t'.join(header)
    for i in range(file_df.shape[0]):
        #print file_df[0][i]
        sample    = re.split(r'_', file_df['Sample'][i])[0]
        timepoint = re.split(r'_', file_df['Sample'][i])[1]
        day       = -1
        if timepoint == 'S':
            day = 0
        elif timepoint == 'M':
            day = 14
        elif timepoint == 'E':
            day = 180
        meta_new = [file_df['Cell.ID'][i], file_df['orig.ident'][i], file_df['nCount_RNA'][i], file_df['nFeature_RNA'][i], sample, timepoint, day, file_df['Celltype'][i], file_df['Celltype_subtype'][i], file_df['Platform'][i]]
        meta_new.extend(clinical_data[file_df['orig.ident'][i]])
        print >> ofile, '\t'.join(map(str, meta_new))
        #print data[file_df['Cell.ID'][i]]
    ofile.close()



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta')
    parser.add_argument('--clinical')
    parser.add_argument('--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.meta) > 0 and len(args.clinical) > 0 
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = "FEL001046_scRNA.metadata.clinical.txt"


    clinical_data = read_clinical(args.clinical) 
    #read cell annot info into a dict
    read_cell_anno(args.meta, clinical_data, args.output)

if __name__ == '__main__':
    main()
