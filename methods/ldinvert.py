from __future__ import print_function, division
import numpy as np
import argparse
import pickle
from estimator import Estimator
from primitives import Dataset, GenomicSubset, SnpSubset
from sparse import ldmatrix


class MLE(Estimator):
    parser = argparse.ArgumentParser()
    parser.add_argument('--Nref', type=int, required=True,
            help='the size of the reference panel to use')
    parser.add_argument('--ld_bandwidth', type=int, required=True,
            help='the maximal ld bandwidth to allow, in SNPs')
    parser.add_argument('--region', type=str, required=True,
            help='the name of the subset of the genome whose heritability should be \
                    analyzed. these files are in data/genome_subsets')

    def readable_name(self):
        return 'MLE,A={},Nref={},ldband={}'.format(
                self.params.region,
                self.params.Nref,
                self.params.ld_bandwidth)

    def preprocessing_folder(self):
        return 'pre.covariance.Nref={}.ldbandwidth={}'.format(
                self.params.Nref,
                self.params.ld_bandwidth)

    def covariance_ldmatrix_file(self, sim, mode='rb'):
        return open(self.path_to_preprocessed_data(sim) + 'covariance.ldmatrix', mode)

    def preprocess(self, sim):
        np.random.seed(0) # to make sure it's the same reference panel every time
        d = Dataset(sim.dataset)
        indivs = d.random_indivs(self.params.Nref)
        ld = ldmatrix.LdMatrix(d, indivs, self.params.ld_bandwidth)
        pickle.dump(ld, self.covariance_ldmatrix_file(sim, mode='wb'), 2)

    def run(self, beta_num, sim):
        R = pickle.load(self.covariance_ldmatrix_file(sim))

        # compute RA from reference panel
        np.random.seed(0) # to make sure it's the same reference panel every time
        d = Dataset(sim.dataset)
        indivs = d.random_indivs(self.params.Nref)
        ss = SnpSubset(GenomicSubset(self.params.region), d.snp_coords())
        RA = ldmatrix.LdMatrix(d, indivs, self.params.ld_bandwidth, snpset_irs=ss.irs)

        # compute the results
        results = []
        for alphahat in sim.sumstats_files(beta_num):
            betahat = R.solve_banded(alphahat)
            results.append(self.compute_statistic(
                alphahat, R, RA, sim.sample_size, self.params.Nref, memoize_bias=True))
            print(results[-1])

        return results

    def compute_statistic(self, alphahat, R, RA, N, Nref, memoize_bias=False):
        if not memoize_bias or not hasattr(self, 'bias'):
            self.bias = R.trace_of_inverse_times_matrix(RA) / N
        betahat = R.solve_banded(alphahat)
        return betahat.dot(RA.covcsr.dot(betahat)) - self.bias

class MLE_reg(MLE):
    parser = argparse.ArgumentParser(add_help=False, parents=[MLE.parser])
    parser.add_argument('--Lambda', type=float, required=True,
            help='the value of lambda by which to regularize LD estimates')

    def readable_name(self):
        return 'MLEreg,A={},Nref={},ldband={},L={}'.format(
                self.params.region,
                self.params.Nref,
                self.params.ld_bandwidth,
                self.params.Lambda)

    def compute_statistic(self, alphahat, R, RA, N, Nref, memoize_bias=False):
        #TODO: should we regularize RA?
        R.add_ridge_and_renormalize(self.params.Lambda)

        if not memoize_bias or not hasattr(self, 'bias'):
            self.bias = R.trace_of_inverse_times_matrix(RA) / N
        betahat = R.solve_banded(alphahat)
        return betahat.dot(RA.covcsr.dot(betahat)) - self.bias

