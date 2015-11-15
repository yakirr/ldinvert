from __future__ import print_function, division
import numpy as np
import pickle
from pysnptools.snpreader import Bed
import sparse.blockdiag as bd
import genome.utils as gutils
import pyutils.iter as it
import hyperparams as hp

hp.load()
pathway_indexset = pickle.load(hp.pathway_file())
#TODO: make the window be in genome coordinates not snps
ldscore_window_in_snps = 500

# sample random set of 14k individuals (this is done for ldscore and for covariance)
np.random.seed(0)
iids = np.random.choice(hp.dataset.N, size=14000, replace=False)

ldscores_A = np.zeros(hp.dataset.M())
ldscores_G = np.zeros(hp.dataset.M())
def compute_ldscores_for_slice(s):
    global ldscores_A, ldscores_G
    print(s)
    genotypes = hp.dataset.genotypes_bedfile()[iids, s[0]:s[1]].read()
    genotypes.standardize(); genotypes.standardize()

    # figure out which snps to compute ld scores for in this slice
    indices = (0 if s[0] == 0 else ldscore_window_in_snps,
        hp.dataset.M() - s[0] if s[1] == hp.dataset.M() else
            genotypes.sid_count-ldscore_window_in_snps)

    one_over_N = 1 / genotypes.iid_count
    def compute_ldscores_for_snp(m):
        start = max(0, m - ldscore_window_in_snps)
        end = min(genotypes.sid_count, m + ldscore_window_in_snps)

        # compute the scores
        window = genotypes.val[:, start:end]
        r2_to_snps_in_window = (genotypes.val[:,m].dot(window) / genotypes.iid_count)**2
        r2_to_snps_in_window -= one_over_N

        # add to appropriate ldscore buckets
        pathway_indicator = np.array([s[0] + m_prime in pathway_indexset
            for m_prime in range(start, end)])
        ldscores_A[s[0] + m] += np.sum(r2_to_snps_in_window * pathway_indicator)
        ldscores_G[s[0] + m] += np.sum(r2_to_snps_in_window * np.logical_not(pathway_indicator))

        # def add_ld_to_neighboring_snp(m_prime):
        #     ldscore = (genotypes.val[:,m].dot(genotypes.val[:,m_prime])
        #             / genotypes.iid_count)**2 - one_over_N
        #     if s[0] + m_prime in pathway_indexset:
        #         ldscores_A[s[0] + m] += ldscore
        #     else:
        #         ldscores_G[s[0] + m] += ldscore
        # map(add_ld_to_neighboring_snp, range(start, end))
    map(compute_ldscores_for_snp, it.show_progress(range(indices[0], indices[1])))
map(compute_ldscores_for_slice, gutils.slices(hp.dataset, buffer_size=ldscore_window_in_snps))

# write output
with hp.ldscores_file_pathway(mode='wb') as ldscores_A_file:
    pickle.dump(ldscores_A, ldscores_A_file, 2)
with hp.ldscores_file_notpathway(mode='wb') as ldscores_G_file:
    pickle.dump(ldscores_G, ldscores_G_file, 2)
