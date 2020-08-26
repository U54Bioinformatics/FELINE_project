#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import re
import os
import argparse
import glob
import psutil


def usage():
    test="name"
    message='''
python scRNA_Seurat_05individual_from_individual_merge_matrix.py --list FEL011027_responder_vs_non_responder.list

--list: list of patients that need to be merged

    '''
    print message

class StringConverter(dict):
    def __contains__(self, item):
        return True

    def __getitem__(self, item):
        return str

    def get(self, default=None):
        return str

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

    print "merging"
    print mem_usage(df)
    #print memory_usage_psutil()
    sys.stdout.flush() 

    # Create new list with merged dataframe
    dfl = []
    dfl.append(df)

    # Join lists
    dfl = dfl + dfs[2:] 
    dfm = mergefiles(dfl, on)
    return dfm


def memory_usage_psutil():
    # return the memory usage in percentage like top
    process = psutil.Process(os.getpid())
    mem = process.get_memory_info()[0] / float(2 ** 20)
    return mem

# We're going to be calculating memory usage a lot,
# so we'll create a function to save us some time!
def mem_usage(pandas_obj):
    if isinstance(pandas_obj,pd.DataFrame):
        usage_b = pandas_obj.memory_usage(deep=True).sum()
    else: # we assume if not a df it's a series
        usage_b = pandas_obj.memory_usage(deep=True)
    usage_mb = usage_b / 1024 ** 2 # convert bytes to megabytes
    return "{:03.2f} MB".format(usage_mb)

#int8	Byte (-128 to 127)
#int16	Integer (-32768 to 32767)
#int32	Integer (-2147483648 to 2147483647)
#int64	Integer (-9223372036854775808 to 9223372036854775807)
#float16	Half precision float: sign bit, 5 bits exponent, 10 bits mantissa
#float32	Single precision float: sign bit, 8 bits exponent, 23 bits mantissa
#float64	Double precision float: sign bit, 11 bits exponent, 52 bits mantissa
def mem_usage_dtype(pandas_obj):
    gl = pandas_obj
    for dtype in ['float','int','object']:
        selected_dtype = gl.select_dtypes(include=[dtype])
        mean_usage_b = selected_dtype.memory_usage(deep=True).sum()
        mean_usage_mb = mean_usage_b / 1024 ** 2
        print("Total memory usage for {} columns: {:03.2f} MB".format(dtype,mean_usage_mb))
        #converted dtype to reduce memory useage
        if dtype == 'float':
            converted_float = selected_dtype.astype('float32')
            #replacing orginal data in df
            #gl.loc[:,selected_dtype.columns] = converted_float
            mean_usage_b = converted_float.memory_usage(deep=True).sum()
            mean_usage_mb = mean_usage_b / 1024 ** 2
            print("Total memory usage for {} columns as float32: {:03.2f} MB".format(dtype,mean_usage_mb))
        if dtype == 'int':
            converted_int = selected_dtype.astype('int16')
            #gl.loc[:,selected_dtype.columns] = converted_int
            mean_usage_b = converted_int.memory_usage(deep=True).sum()
            mean_usage_mb = mean_usage_b / 1024 ** 2
            print("Total memory usage for {} columns as int16: {:03.2f} MB".format(dtype,mean_usage_mb))
    return gl

def read_file_pd(infile):
    # determine and optimize dtype
    # Sample 100 rows of data to determine dtypes.
    file_test = pd.read_csv(infile, sep="\t", header=0, nrows=100)
    float_cols = [c for c in file_test if file_test[c].dtype == "float64"]
    int_cols = [c for c in file_test if file_test[c].dtype == "int64"]
    if float_cols > 0:
        dtype_cols = {c: np.float16 for c in float_cols}
    elif int_cols > 0:
        dtype_cols = {c: np.int16 for c in int_cols}
    file_df = pd.read_csv(infile, sep="\t", header=0, dtype=dtype_cols)
    # check memory usage
    print "infile: %s" %(infile)
    print "original size"
    print file_df.dtypes[1:5,]
    print file_df.iloc[1:5,1:50]
    print mem_usage(file_df)
    #file_df_mini = mem_usage_dtype(file_df)
    #print "converted size"
    #print file_df_mini.dtypes[1:5,]
    #print file_df_mini.iloc[1:5,1:50]
    #print mem_usage(file_df_mini)
    sys.stdout.flush()
    #return file_df_mini
    return file_df

def get_file_df(patients, filetype):
    df_list   = []
    for p in patients:
        filename = '%s%s' %(p, filetype)
        if os.path.exists(filename):
            df_file = read_file_pd(filename)
            df_list.append(df_file)
    return df_list

#FEL027_M_TATATCCTCAAGGTGG
def read_patient_list(infile, prefix):
    file_df  = pd.read_table(infile, header=None, names=['Patient'])
    patients = ['%s_%s' %(prefix, x) for x in list(file_df['Patient'])]
    print patients
    return patients


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--list')
    parser.add_argument('--prefix')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.list) > 0 and len(args.prefix) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = "test_merge"

    #file_list = glob.glob('%s/*.ssGSEA.scores.txt' %(args.input))
    #df_list   = []
    #for f in sorted(file_list):
    #    df_file = read_file_pd(f)
    #    df_list.append(df_file)

    patients = read_patient_list(args.list, args.prefix)

    #filetypes = ['_10x.seurat.cell_type.anno.txt']
    #for ft in filetypes:
    #    outfile = '%s.%s' %(args.list, ft)
    #    if re.compile(r'.list$').search(args.list):
    #        outfile = re.sub(r'.list', r'%s' %(ft), args.list)
    #    df_list = get_file_df(patients, ft)
    #    df_merged = pd.concat(df_list)
    #    df_merged.to_csv(outfile, sep="\t", index=False)


    #filetypes = ['_10x_gene_symbols.scaled.counts.normalized_to_control.txt', '_10x_gene_symbols.scaled.Cancer_cells.mean.txt', '_10x_gene_symbols.scaled.Endothelial_cells.mean.txt', '_10x_gene_symbols.scaled.Fibroblasts.mean.txt', '_10x_gene_symbols.scaled.Macrophages.mean.txt']
    #filetypes = ['_gene_symbols.normalized.counts.txt', '_gene_symbols.scaled.counts.txt', '_gene_symbols.CPM.txt', '_gene_symbols.raw.counts.txt']
    #filetypes = ['_gene_symbols.normalized.counts.txt']
    filetypes = ['_gene_symbols.raw.counts.txt', '_gene_symbols.CPM.txt']
    #filetypes = ['_gene_symbols.raw.counts.txt']
    #filetypes = ['_gene_symbols.CPM.txt']
    for ft in filetypes:
        print ft
        outfile = '%s%s' %(args.prefix, ft)
        df_list = get_file_df(patients, ft)
        df_merged = mergefiles(df_list, on="Gene.ID") 
        df_merged.to_csv(outfile, sep="\t", index=False)
    #merge
    #outfile   = '%s.ssGSEA.scores.txt' %(args.output)
    #df_merged = mergefiles(df_list, on="Gene Set")
    #df_merged.to_csv(outfile, sep="\t", index=False)

if __name__ == '__main__':
    main()

