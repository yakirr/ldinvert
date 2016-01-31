from __future__ import print_function, division
import numpy as np
import argparse
import pickle
from estimator import Estimator
from primitives import Dataset, GenomicSubset, SnpSubset
# from sparse import ldmatrix
from sparse.blockdiag import BlockDiag


class Truth(Estimator):
    parser = argparse.ArgumentParser()
    parser.add_argument('--ld_bandwidth', type=float, required=True,
            help='the maximal ld bandwidth to allow, in Morgans')
    parser.add_argument('--region', type=str, required=True,
            help='the name of the subset of the genome whose heritability should be \
                    analyzed. these files are in data/genome_subsets')

    def readable_name(self):
        return 'Truth,A={},ld_band={}'.format(
                self.params.region,
                self.params.ld_bandwidth)

    def preprocessing_folder(self):
        return 'pre.covariance.A={}.Nref=full.ldbandwidth={}'.format(
                self.params.region,
                self.params.ld_bandwidth)

    def RA_file(self, sim, mode='rb'):
        return open(self.path_to_preprocessed_data(sim) + 'RA.bd', mode)

    def preprocess(self, sim):
        d = Dataset(sim.dataset)
        gs = GenomicSubset(self.params.region)
        ss = SnpSubset(gs, d)
        RA = BlockDiag.ld_matrix(d, ss.irs, self.params.ld_bandwidth)
        pickle.dump(RA, self.RA_file(sim, mode='wb'), 2)

    def run(self, beta_num, sim):
        RA = pickle.load(self.RA_file(sim))
        beta = pickle.load(sim.beta_file(beta_num))
        beta = BlockDiag.from_big1darray(beta, RA.irs)
        results = [beta.dot(RA.dot(beta))]
        print(results[-1])
        return results
