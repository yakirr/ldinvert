from __future__ import print_function, division
import numpy as np
import argparse
import pickle
from estimator import Estimator
from primitives import Dataset, GenomicSubset, SnpSubset
# from sparse import ldmatrix
from sparse.blockdiag import BlockDiag


class MLE(Estimator):
    parser = argparse.ArgumentParser()
    parser.add_argument('--Nref', type=int, required=True,
            help='the size of the reference panel to use')
    parser.add_argument('--ld_bandwidth', type=float, required=True,
            help='the maximal ld bandwidth to allow, in Morgans')
    parser.add_argument('--region', type=str, required=True,
            help='the name of the subset of the genome whose heritability should be \
                    analyzed. these files are in data/genome_subsets')

    def readable_name(self):
        return 'MLE,A={},Nref={},ldband={}'.format(
                self.params.region,
                self.params.Nref,
                self.params.ld_bandwidth)

    def preprocessing_folder(self):
        return 'pre.covariance.A={}.Nref={}.ldbandwidth={}'.format(
                self.params.region,
                self.params.Nref,
                self.params.ld_bandwidth)

    def R_file(self, sim, mode='rb'):
        return open(self.path_to_preprocessed_data(sim) + 'R.bd', mode)

    def RA_file(self, sim, mode='rb'):
        return open(self.path_to_preprocessed_data(sim) + 'RA.bd', mode)

    def preprocess(self, sim):
        np.random.seed(0) # to make sure it's the same reference panel every time
        d = Dataset(sim.dataset)
        indivs = d.random_indivs(self.params.Nref)
        gs = GenomicSubset(self.params.region)
        ss = SnpSubset(gs, d)
        buffered_ss = ss.expanded_by(self.params.ld_bandwidth)
        R = BlockDiag.ld_matrix(d, buffered_ss.irs, self.params.ld_bandwidth, indivs=indivs)
        RA = R.zero_outside_irs(ss.irs)
        pickle.dump(R, self.R_file(sim, mode='wb'), 2)
        pickle.dump(RA, self.RA_file(sim, mode='wb'), 2)

    def run(self, beta_num, sim):
        R = pickle.load(self.R_file(sim))
        RA = pickle.load(self.RA_file(sim))

        # compute the results
        results = []
        for alphahat in sim.sumstats_files(beta_num):
            alphahat = BlockDiag.from_big1darray(alphahat, R.irs)
            results.append(self.compute_statistic(
                alphahat, R, RA, sim.sample_size, self.params.Nref, memoize=True))
            print(len(results), results[-1])

        return results

    def compute_statistic(self, alphahat, R, RA, N, Nref, memoize=False):
        if not memoize or not hasattr(self, 'bias'):
            print('computing bias')
            self.bias = BlockDiag.solve(R, RA).trace() / N
            print('bias =', self.bias)
        betahat = BlockDiag.solve(R, alphahat)
        return betahat.dot(RA.dot(betahat)) - self.bias

    # def run(self, beta_num, sim):
    #     # compute RA from reference panel
    #     np.random.seed(0) # to make sure it's the same reference panel every time
    #     d = Dataset(sim.dataset)
    #     indivs = d.random_indivs(self.params.Nref)
    #     gs = GenomicSubset(self.params.region)
    #     ss = SnpSubset(gs, d)
    #     RA = ldmatrix.LdMatrix(d, indivs, self.params.ld_bandwidth, snpset_irs=ss.irs)
    #     # R = pickle.load(self.covariance_ldmatrix_file(sim))
    #     gs.expand_by(1)
    #     ss_expanded = SnpSubset(gs, d)
    #     R = ldmatrix.LdMatrix(d, indivs, self.params.ld_bandwidth, snpset_irs=ss_expanded.irs)

    #     # compute the results
    #     results = []
    #     for alphahat in sim.sumstats_files(beta_num):
    #         results.append(self.compute_statistic(
    #             alphahat, R, RA, sim.sample_size, self.params.Nref, memoize=True))
    #         print(len(results), results[-1])

    #     return results

    # def compute_statistic(self, alphahat, R, RA, N, Nref, memoize=False):
    #     if not memoize or not hasattr(self, 'bias'):
    #         print('computing bias')
    #         self.bias = R.trace_of_inverse_times_matrix(RA) / N
    #         print('bias =', self.bias)
    #     betahat = R.solve_banded(alphahat)
    #     return betahat.dot(RA.covcsr.dot(betahat)) - self.bias

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

    def compute_statistic(self, alphahat, R, RA, N, Nref, memoize=False):
        #TODO: should we regularize RA?
        print('regularizing R...')
        Rreg = R.add_ridge(self.params.Lambda, renormalize=True)
        if not memoize or not hasattr(self, 'bias'):
            print('done.computing bias...')
            self.bias = BlockDiag.solve(Rreg, RA).trace() / N
            print('bias =', self.bias)
        betahat = BlockDiag.solve(Rreg, alphahat)
        return betahat.dot(RA.dot(betahat)) - self.bias

    # def compute_statistic(self, alphahat, R, RA, N, Nref, memoize=False):
    #     #TODO: should we regularize RA?
    #     if not memoize or not hasattr(self, 'bias'):
    #         print('regularizing R...')
    #         R.add_ridge_and_renormalize(self.params.Lambda)
    #         print('done.computing bias...')
    #         self.bias = R.trace_of_inverse_times_matrix(RA) / N
    #         print('bias =', self.bias)
    #     betahat = R.solve_banded(alphahat)
    #     return betahat.dot(RA.covcsr.dot(betahat)) - self.bias

