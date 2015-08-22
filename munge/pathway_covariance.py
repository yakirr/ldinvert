from __future__ import print_function, division
import numpy as np
import pickle
from pysnptools.snpreader import Bed
from pybedtools import BedTool, Interval
import sparse.blockdiag as bd
import hyperparams as hp

regions_to_indexsets = pickle.load(hp.pathway_with_flanks_file())

covariance_matrices = {}
for region, indexset in regions_to_indexsets.items():
    print('region', region, end='')

    chrom_data_on_disk = hp.genotypes_bed_file(region.chrom[3:])
    region_data = chrom_data_on_disk[:, indexset].read()
    print('read in all snps in region. shape =', region_data.val.shape)

    print('standardizing...')
    # we need the second call because the first one standardized the variance before substituting 0s for the nans.
    region_data.standardize(); region_data.standardize()

    print('computing covariance...')
    covariance_matrices[region] = region_data.val.T.dot(region_data.val) / region_data.iid_count

    # TODO: regularize the estimated covariance matrix by banding it with a 1MB (or 1cM) window

print('assembling into BlockDiagArray...')
bda = bd.BlockDiagArray(covariance_matrices, regions_to_indexsets)

print('saving...')
pickle.dump(bda, hp.covariance_around_pathway_file(mode='wb'), 2)
