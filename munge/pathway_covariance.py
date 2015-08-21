from __future__ import print_function, division
import numpy as np
import pickle
from pysnptools.snpreader import Bed
from pybedtools import BedTool, Interval
import sparse.blockdiag as bd
import partition.paths as pathto

dataset_name = 'WT1QC_fewnans'
pathway_name = '99'

regions_to_indexsets = pickle.load(
        open(pathto.pathways_with_flanks + pathway_name + '.' + dataset_name + \
                '/merged.regions_to_indexsets', 'rb'))

covariance_matrices = {}
for region, indexset in regions_to_indexsets.items():
    print('region', region, end='')

    chrom_data_on_disk = Bed(pathto.genotypes + dataset_name + '/all.' + region.chrom[3:])
    region_data = chrom_data_on_disk[:, indexset].read()
    print('read in all snps in region. shape =', region_data.val.shape)

    print('standardizing...')
    # we need the second call because the first one standardized the variance before substituting 0s for the nans.
    region_data.standardize()
    region_data.standardize()

    print('computing covariance...')
    covariance_matrices[region] = np.cov(region_data.val.T)

    # TODO: regularize the estimated covariance matrix by banding it with a 1MB (or 1cM) window

print('assembling into BlockDiagArray...')
bda = bd.BlockDiagArray(covariance_matrices, regions_to_indexsets)

print('saving...')
pickle.dump(bda, open(pathto.pathways_with_flanks + pathway_name + '.' + dataset_name + '/merged.covariance.bda', 'wb'), 2)

# for using savez instead of pickling a BlockDiagArray
# np.savez(pathto.pathways_with_flanks + pathway_name + '.' + dataset_name + '/merged.covariance', **covariance_matrices)
