from __future__ import print_function, division
import numpy as np
import pickle
from pysnptools.snpreader import Bed
import sparse.blockdiag as bd
import genome.utils as gutils
import pyutils.iter as it
import hyperparams as hp

hp.load()
pathway = pickle.load(hp.pathway_file())
#TODO: make the window be in genome coordinates not snps
ldscore_window_in_snps = 100

chrnum_to_ldscores_A = {}
chrnum_to_ldscores_G = {}
def compute_ldscores_for_chrom(chrnum):
    global chrnum_to_ldscores_A, chrnum_to_ldscores_G
    print(chrnum)
    genotypes = hp.genotypes_bed_file(chrnum).read()
    genotypes.standardize(); genotypes.standardize()
    chrnum_to_ldscores_A[chrnum] = np.zeros(genotypes.sid_count)
    chrnum_to_ldscores_G[chrnum] = np.zeros(genotypes.sid_count)
    pathway_indexset = gutils.indexset_wrt_chromosome(chrnum, pathway)

    one_over_N = 1 / genotypes.iid_count
    def compute_ldscores_for_snp(m):
        start = max(0, m - int(ldscore_window_in_snps / 2))
        end = min(genotypes.sid_count, m + int(ldscore_window_in_snps / 2))

        def add_ld_to_neighboring_snp(m_prime):
            ldscore = (genotypes.val[:,m].dot(genotypes.val[:,m_prime]) / genotypes.iid_count)**2 - one_over_N
            if m_prime in pathway_indexset:
                chrnum_to_ldscores_A[chrnum][m] += ldscore
            else:
                chrnum_to_ldscores_G[chrnum][m] += ldscore
        map(add_ld_to_neighboring_snp, range(start, end))
    map(compute_ldscores_for_snp, it.show_progress(range(genotypes.sid_count)))
map(compute_ldscores_for_chrom, hp.chromosomes())

# write output
with ldscores_A_file, ldscores_G_file = hp.ldscores_files(mode='wb'):
    pickle.dump(chrnum_to_ldscores_A, ldscores_A_file, 2)
    pickle.dump(chrnum_to_ldscores_G, ldscores_G_file, 2)
