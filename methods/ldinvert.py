from __future__ import print_function, division
import numpy as np
import argparse
import pickle
from estimator import Estimator
from primitives import Dataset, GenomicSubset, SnpSubset
# from sparse import ldmatrix
from sparse.blockdiag import BlockDiag


class MLE(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--ld_bandwidth', type=float, required=True,
            help='the maximal ld bandwidth to allow, in Morgans')
    parser.add_argument('--region', type=str, required=True,
            help='the name of the subset of the genome whose heritability should be \
                    analyzed. these files are in data/genome_subsets')

    def readable_name(self):
        return 'MLE,A={},ref={},ldband={}'.format(
                self.params.region,
                self.params.refpanel,
                self.params.ld_bandwidth)

    def preprocessing_foldername(self):
        return 'pre.covariance.A={}.ldbandwidth={}'.format(
                self.params.region,
                self.params.ld_bandwidth)

    def R_file(self, mode='rb'):
        return open(self.path_to_preprocessed_data() + 'R.bd', mode)

    def RA_file(self, mode='rb'):
        return open(self.path_to_preprocessed_data() + 'RA.bd', mode)

    def preprocess(self):
        gs = GenomicSubset(self.params.region)
        ss = SnpSubset(self.refpanel, bedtool=gs.bedtool)
        buffered_ss = ss.expanded_by(self.params.ld_bandwidth)
        R = BlockDiag.ld_matrix(self.refpanel, buffered_ss.irs, self.params.ld_bandwidth)
        RA = R.zero_outside_irs(ss.irs)
        pickle.dump(R, self.R_file(mode='wb'), 2)
        pickle.dump(RA, self.RA_file(mode='wb'), 2)

    def run(self, beta_num, sim):
        R = pickle.load(self.R_file())
        RA = pickle.load(self.RA_file())

        # compute the results
        results = []
        for alphahat in sim.sumstats_files(beta_num):
            alphahat = BlockDiag.from_big1darray(alphahat, R.irs)
            results.append(self.compute_statistic(
                alphahat, R, RA, sim.sample_size, self.refpanel.N, memoize=True))
            print(len(results), results[-1])

        return results

    def compute_statistic(self, alphahat, R, RA, N, Nref, memoize=False):
        try:
            if not memoize or not hasattr(self, 'bias'):
                print('computing bias')
                self.bias = BlockDiag.solve(R, RA).trace() / N
                print('bias =', self.bias)
            betahat = BlockDiag.solve(R, alphahat)
            return betahat.dot(RA.dot(betahat)) - self.bias
        except np.linalg.linalg.LinAlgError:
            print('R was singular. Its shape was', R.shape(), 'and Nref=', Nref)
            return 0


class MLE_reg(MLE):
    parser = argparse.ArgumentParser(add_help=False, parents=[MLE.parser])
    parser.add_argument('--Lambda', type=float, required=True,
            help='the value of lambda by which to regularize LD estimates')

    def readable_name(self):
        return 'MLEreg,A={},ref={},ldband={},L={}'.format(
                self.params.region,
                self.params.refpanel,
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


# version of MLE_reg that doesn't band when estimating LD around a category.
class MLE_reg_noband(MLE_reg):
    def readable_name(self):
        return 'MLE_reg_noband,A={},ref={},ldband={}'.format(
                self.params.region,
                self.params.refpanel,
                self.params.ld_bandwidth)

    def preprocessing_foldername(self):
        return 'pre.unbanded_covariance.A={}.ldbandwidth={}'.format(
                self.params.region,
                self.params.ld_bandwidth)

    def preprocess(self):
        gs = GenomicSubset(self.params.region)
        ss = SnpSubset(self.refpanel, bedtool=gs.bedtool)
        buffered_ss = ss.expanded_by(self.params.ld_bandwidth)
        R = BlockDiag.ld_matrix(self.refpanel, buffered_ss.irs, 1000000) # bandwidth=infty
        RA = R.zero_outside_irs(ss.irs)
        pickle.dump(R, self.R_file(mode='wb'), 2)
        pickle.dump(RA, self.RA_file(mode='wb'), 2)


if __name__ == '__main__':
    import primitives.genome, primitives.simulation, primitives.dataset
    reload(primitives.genome)
    reload(primitives.simulation)
    reload(primitives.dataset)
    from primitives.genome import GenomicSubset, SnpSubset
    from primitives.simulation import SumstatSimulation
    from primitives.dataset import Dataset
    est = MLE(refpanel='tinyGERA_ref14k', region='tiny', ld_bandwidth=0.01)
    # est = MLE_reg_noband(refpanel='tinyGERA_ref2k', region='tiny', ld_bandwidth=0.01,
            # Lambda=0.01)
    sim = SumstatSimulation('tinyGERA.tiny_inf')
    # est.preprocess()
    est.run(1, sim)