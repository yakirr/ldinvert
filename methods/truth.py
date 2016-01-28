from __future__ import print_function, division
import numpy as np
import argparse
import pickle
from estimator import Estimator
from primitives import Dataset, GenomicSubset, SnpSubset
from sparse import ldmatrix


class Truth(Estimator):
    parser = argparse.ArgumentParser()
    parser.add_argument('--ld_bandwidth', type=int, required=True,
            help='the maximal ld bandwidth to allow, in SNPs')
    parser.add_argument('--region', type=str, required=True,
            help='the name of the subset of the genome whose heritability should be \
                    analyzed. these files are in data/genome_subsets')

    def readable_name(self):
        return 'Truth,A={},ld_band={}'.format(
                self.params.region,
                self.params.ld_bandwidth)

    def preprocessing_folder(self):
        return 'empty'

    def preprocess(self, sim):
        pass

    def run(self, beta_num, sim):
        # compute RA
        d = Dataset(sim.dataset)
        indivs = np.arange(d.N)
        ss = SnpSubset(GenomicSubset(self.params.region), d.snp_coords())
        RA = ldmatrix.LdMatrix(d, indivs, self.params.ld_bandwidth, snpset_irs=ss.irs)

        # compute the results
        results = []
        beta = pickle.load(sim.beta_file(beta_num))
        results.append(beta.dot(RA.covcsr.dot(beta)))
        print(results[-1])

        return results
