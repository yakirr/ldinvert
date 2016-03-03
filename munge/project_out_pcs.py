from __future__ import print_function, division
import pandas as pd
import numpy as np
import argparse
from pysnptools.snpreader import Bed, SnpData

parser = argparse.ArgumentParser()
parser.add_argument('--bfile', type=str, required=True,
        help='the bfile to project the PCs out of')
parser.add_argument('--pcs', type=str, required=True,
        help='path to the pcs file')

args = parser.parse_args()
pcs_df = pd.read_csv(args.pcs, sep='\t', usecols=[0,1,5,6,7,8,9,10,11,12,13,14],
        # names=['FID','IID','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'],
        header=0)

fam_df = pd.read_csv(args.bfile + '.fam', delim_whitespace=True,
        usecols=[0,1], names=['FID','IID'])

merged_df = pd.merge(fam_df, pcs_df, on='IID').set_index('IID')
pcs = merged_df.ix[fam_df.ix[:,1], 2:]
Q, R = np.linalg.qr(pcs)

print('reading dataset')
dataset = Bed(args.bfile).read().standardize()
dataset.standardize()

import pdb; pdb.set_trace()
Bed.write('temp', dataset)

print('projecting data')
X_Q = Q.T.dot(dataset.val)

print('unprojecting')
X_rr = Q.dot(X_Q)

print('subtracting out population structure')
X = dataset.val - X_rr

print('writing')
newbed = SnpData(dataset.iid, dataset.sid, X, pos=dataset.pos)
Bed.write('temp.bed', newbed)
