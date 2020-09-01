import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import os
import sys

prefix=''

if len(sys.argv) == 1:
    print('input file prefix: python scrublet_doublet.py FEL011_S')
    exit()
else:
    prefix=sys.argv[1]
    print(prefix)

#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42


counts_matrix = scipy.io.mmread(prefix + '.matrix.mtx').T.tocsc()
genes = np.array(scr.load_genes(prefix + '.genes.tsv', delimiter='\t', column=0))
cells = pd.read_table(prefix + '.barcodes.tsv', header=None)
cells.columns = ["Cell.ID"]



print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

#indexnames  = list(counts_matrix.index)
#columnnames = list(counts_matrix.columns)

#print('10 index values: {}'.format(indexnames[1:10]))
#print('10 column values: {}'.format(counts_matrix[1:3,1:3]))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)


doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=1, 
                                                          min_cells=100, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

predicted_doublets_025 = scrub.call_doublets(threshold=0.25)

out_df = pd.DataFrame({
    'Cell.ID' : cells['Cell.ID'],
    'scrublet_doublet_score' : doublet_scores,
    'scrublet_doublet_call1' : predicted_doublets,
    'scrublet_doublet_call2' : predicted_doublets_025
    })

out_df.to_csv('%s.scrublet_out.csv' %(prefix), sep="\t", index=False)

with PdfPages('%s.scrublet_out.pdf' %(prefix)) as pdf:
    scrub.plot_histogram();
    pdf.savefig()


    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    # # Uncomment to run tSNE - slow
    # print('Running tSNE...')
    # scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))
    # # Uncomment to run force layout - slow
    # print('Running ForceAtlas2...')
    # scrub.set_embedding('FA', scr.get_force_layout(scrub.manifold_obs_, n_neighbors=5. n_iter=1000))
    pdf.savefig()
    print('Done.')

    scrub.plot_embedding('UMAP', order_points=True);
    # scrub.plot_embedding('tSNE', order_points=True);
    # scrub.plot_embedding('FA', order_points=True);
    pdf.savefig()

