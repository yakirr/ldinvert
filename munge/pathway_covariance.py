from __future__ import print_function, division
import numpy as np
import pickle
from pysnptools.snpreader import Bed
from pybedtools import BedTool, Interval
import sparse.blockdiag as bd
import genome.utils as gutils
import hyperparams as hp

hp.load()
merged_intrangeset = pickle.load(hp.pathway_with_flanks_file())

# sample random set of 14k individuals (this is done for ldscore and for covariance)
np.random.seed(0)
iids = np.random.choice(hp.dataset.N, size=14000, replace=False)

covariance_matrices = {}
for r in merged_intrangeset.ranges():
    print('range', r)

    range_data_on_disk = hp.dataset.genotypes_bedfile()
    range_data = range_data_on_disk[iids, r[0]:r[1]].read()
    print('read in all snps in region. shape =', range_data.val.shape)

    print('standardizing...')
    # we need the second call because the first one standardized the variance
    # before substituting 0s for the nans.
    range_data.standardize(); range_data.standardize()

    print('computing covariance...')
    covariance_matrices[r] = \
            range_data.val.T.dot(range_data.val) / range_data.iid_count

    # TODO: regularize estimated covariance matrix by banding it with a 1MB (or 1cM) window

print('assembling into BlockDiagArray...')
bda = bd.BlockDiagArray(covariance_matrices)

print('saving...')
pickle.dump(bda, hp.covariance_around_pathway_file(mode='wb'), 2)
