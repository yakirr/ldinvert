from __future__ import print_function
import pandas as pd
import numpy as np
import gzip

# sumstats = pd.read_csv('/groups/price/yakir/data/sumstats/CD.sumstats.gz',
#         header=0, sep='\t', compression='gzip')
for chrnum in range(1,23):
    print(chrnum)
    annot = pd.read_csv(
            '/groups/price/yakir/data/signed_annot/T_cell_Andersson_Enhancer.' + \
                str(chrnum) + '.annot.gz',
            header=0, sep='\t', compression='gzip')
    annotname = 'T_cell_Andersson_Enhancer'

    merged = annot.merge(sumstats, how='left', on=['SNP'])
    merged.ix[merged['Z'] < 0, annotname] *= -1
    annot = merged[annot.columns.values]
    print('writing')
    with gzip.open('/groups/price/yakir/data/signed_annot/T_cell_Andersson_Enhancer.cheat.'+\
                str(chrnum) + '.annot.gz', 'wt') as f:
        annot.to_csv(f, index=False, sep='\t')
