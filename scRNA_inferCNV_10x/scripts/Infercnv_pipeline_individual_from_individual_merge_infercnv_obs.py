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
python Infercnv_pipeline_individual_from_individual_merge_infercnv_obs.py --subsample FEL028P131a,FEL028P131b

Merge infercnv results from FEL028P131a, FEL028P131b into FEL028P131.
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
    file_df = pd.read_csv(infile, sep=" ", header=0)
    file_df["Gene.ID"] = list(file_df.index.values)
    return file_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subsample')
    parser.add_argument('--samplesize')
    parser.add_argument('--cells')
    parser.add_argument('--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.subsample) > 0
    except:
        usage()
        sys.exit(2)

    if not args.samplesize:
        args.samplesize = "all"

    if not args.output:
        args.output = 'infercnv_merge'

    #parse and prepare input and output folders
    patient     = ''
    subsets     = re.split(',', args.subsample)
    subsets_out = []
    suffix = 'output_dir_plot'
    if args.cells:
        suffix = '%s_%scells' %(suffix, args.cells)
    for subset in subsets:
        out = '%s_10x_counts.%s.%s' %(subset, args.samplesize, suffix) 
        subsets_out.append(out)
        args.output = '%s_10x_counts.%s.%s' %(subset[:-1], args.samplesize, suffix)
        patient     = subset[:-1]     

    #make output folder
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    print 'Patient: %s' %(patient)
    print 'Input folders: %s' %(subsets_out)
    print 'Output folders: %s' %(args.output)

    #merge observation files
    obsfiles = ['infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt', 'infercnv.observations.txt', 'infercnv.preliminary.observations.txt']
    #obsfiles = ['test.data.txt']
    for subsets_obs in obsfiles:
        subsets_obs_df_list = []
        for subsets_folder in subsets_out:
            file_name = '%s/%s' %(subsets_folder, subsets_obs)
            file_df   = read_file_pd(file_name)
            subsets_obs_df_list.append(file_df)
        subsets_obs_df_merged = mergefiles(subsets_obs_df_list, on="Gene.ID")
        subsets_obs_df_merged.set_index("Gene.ID", inplace=True, drop=True)
        #write output into file with orignial name
        subsets_obs_outfile1  = '%s/%s' %(args.output, subsets_obs) 
        subsets_obs_df_merged.to_csv(subsets_obs_outfile1, sep=" ", index=True, index_label=False)
        #write output into file with patient id as prefix
        subsets_obs_outfile2  = '%s/%s.%s' %(args.output, patient, subsets_obs)
        subsets_obs_df_merged.to_csv(subsets_obs_outfile2, sep=" ", index=True, index_label=False)

    #file_list = glob.glob('%s/*.ssGSEA.scores.txt' %(args.input))
    #df_list   = []
    #for f in sorted(file_list):
    #    df_file = read_file_pd(f)
    #    df_list.append(df_file)
    
    #merge
    #outfile   = '%s.ssGSEA.scores.txt' %(args.output)
    #df_merged = mergefiles(df_list, on="Gene Set")
    #df_merged.to_csv(outfile, sep="\t", index=False)

if __name__ == '__main__':
    main()

