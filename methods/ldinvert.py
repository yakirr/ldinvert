from __future__ import print_function, division
import numpy as np
import matplotlib
import argparse
import pickle
from estimator import Estimator
from primitives import Dataset, GenomicSubset, SnpSubset
# from sparse import ldmatrix
from sparse.blockdiag import BlockDiag


class MLE(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--ld_window', type=int, required=True,
            help='the maximal ld window to allow, in mM')
    parser.add_argument('--region', type=str, required=True,
            help='the name of the subset of the genome whose heritability should be \
                    analyzed. these files are in data/genome_subsets')

    def readable_name(self):
        return 'MLE,A={},ref={},ldwindow={}'.format(
                self.params.region,
                self.params.refpanel,
                self.params.ld_window)

    def preprocessing_foldername(self):
        return 'pre.covariance.A={}.ldwindow={}'.format(
                self.params.region,
                self.params.ld_window)

    def R_file(self, mode='rb'):
        return open(self.path_to_preprocessed_data() + 'R.bd', mode)
    def R_plotfilename(self):
        return self.path_to_preprocessed_data() + 'R.png'

    def RA_file(self, mode='rb'):
        return open(self.path_to_preprocessed_data() + 'RA.bd', mode)

    def preprocess(self):
        matplotlib.use('Agg')
        gs = GenomicSubset(self.params.region)
        A = SnpSubset(self.refpanel, bedtool=gs.bedtool)
        W = A.expanded_by(self.params.ld_window / 1000.)
        R = BlockDiag.ld_matrix(self.refpanel, W.irs.ranges(), 300, band_units='SNPs')
        pickle.dump(R, self.R_file(mode='wb'), 2)
        # R.plot(A.irs, filename=self.R_plotfilename())
        RA = R.zero_outside_irs(A.irs)
        pickle.dump(RA, self.RA_file(mode='wb'), 2)

    def run(self, beta_num, sim):
        R = pickle.load(self.R_file())
        RA = pickle.load(self.RA_file())

        # compute the results
        results = []
        for alphahat in sim.sumstats_aligned_to_refpanel(beta_num, self.refpanel):
            alphahat = BlockDiag.from_big1darray(alphahat, R.ranges())
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
        return 'MLEreg,A={},ref={},ldwindow={},L={}'.format(
                self.params.region,
                self.params.refpanel,
                self.params.ld_window,
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
class MLE_rnb(MLE_reg):
    parser = argparse.ArgumentParser(add_help=False, parents=[MLE_reg.parser])
    parser.add_argument('--units', type=str, required=True,
            help='either mM or bp. Sets how the estimator chooses to draw the LD \
                    window around the category')
    parser.add_argument('--Radjust', type=str, required=False, default='none',
            help='whether to bias adjust the estimate of Rinv for finite reference panel \
                    sample size. "before" means adjust before regularization \
                    and inversion. "after" means adjust after regularization and \
                    inversion. Anything else means don\'t adjust')
    parser.add_argument('--RAreg', type=int, required=False, default=0,
            help='whether to regularize RA in addition to Rinv. 0 means no, 1 means yes.')
    parser.add_argument('--sigma2g', type=float, required=False, default=0,
            help='the value to use for the bias-correction corresponding to a finite pop \
                    size')
    parser.add_argument('--pop_size', type=int, required=False, default=54734,
            help='the size of the population that individuals are sampled from')
    parser.add_argument('--avgunbiased', type=int, required=False, default=0,
            help='1 to apply the average bias correction using the full GERA \
                    genotypes, 0 otherwise.')
    parser.add_argument('--prune_regions', type=int, required=False, default=0,
            help='the number of regions to prune in order to reduce variance')

    def readable_name(self):
        result = 'MLE_rnb,A={},ref={},ldwindow={},L={},units={},Radjust={}'.format(
                    self.params.region,
                    self.params.refpanel,
                    self.params.ld_window,
                    self.params.Lambda,
                    self.params.units,
                    self.params.Radjust) + \
                'RAreg={},sigma2g={:0.2f},pop_size={}'.format(
                    self.params.RAreg,
                    self.params.sigma2g,
                    self.params.pop_size)
        if self.params.avgunbiased:
            result += 'avgunbiased'
        if self.params.prune_regions > 0:
            result +=',prune={}'.format(self.params.prune_regions)
        return result

    def window(self, A):
        if self.params.units == 'mM':
            units = 'Morgans'
            window = self.params.ld_window / 1000.
        elif self.params.units == 'bp':
            units = 'bp'
            window = self.params.ld_window
        else:
            print('ERROR: units must be either mM or bp')
        return A.expanded_by(window, units=units)

    def preprocessing_foldername(self):
        return 'pre.unbanded_covariance.A={}.ldwindow={}.units={}'.format(
                self.params.region,
                self.params.ld_window,
                self.params.units)

    def preprocess(self):
        matplotlib.use('Agg')
        gs = GenomicSubset(self.params.region)
        A = SnpSubset(self.refpanel, bedtool=gs.bedtool)
        W = self.window(A)
        R = BlockDiag.ld_matrix(self.refpanel, W.irs.ranges(), 1000000) # bandwidth=infty
        pickle.dump(R, self.R_file(mode='wb'), 2)
        try: # if the plotting has some error we don't want to not save the stuff
            # R.plot(A.irs, filename=self.R_plotfilename())
            pass
        except:
            pass
        RA = R.zero_outside_irs(A.irs)
        pickle.dump(RA, self.RA_file(mode='wb'), 2)

    def run(self, beta_num, sim):
        R = pickle.load(self.R_file())
        RA = pickle.load(self.RA_file())

        if self.params.prune_regions > 0:
            def var(L, LA, h2A, N):
                LinvLA = np.linalg.solve(L, LA)
                tr1 = np.einsum('ij,ji', LinvLA, LinvLA)
                tr2 = np.einsum('ij,ji', LA, LinvLA)
                return 2*tr1/float(N)**2 + 2*tr2*h2A/(float(N) * float(750))

            print('computing variances')
            variances = {}
            for r in R.ranges():
                variances[r] = var(R.ranges_to_arrays[r], RA.ranges_to_arrays[r],
                        0.05, sim.sample_size)
            print('total variance:', sum(variances.values()))

            sortedrs = R.ranges()
            sortedrs.sort(key=lambda r:variances[r])
            worstrs = sortedrs[-self.params.prune_regions:]
            for r in worstrs:
                print('removing', r)
                del R.ranges_to_arrays[r]
                del RA.ranges_to_arrays[r]
            print('new variance:', sum([variances[r] for r in R.ranges()]))

        print(len(R.ranges()))
        print(len(RA.ranges()))
        # compute the results
        results = []
        for alphahat in sim.sumstats_aligned_to_refpanel(beta_num, self.refpanel):
            alphahat = BlockDiag.from_big1darray(alphahat, R.ranges())
            results.append(self.compute_statistic(
                alphahat, R, RA, sim.sample_size, self.refpanel.N, memoize=True))
            print(len(results), results[-1])

        return results

    def compute_statistic(self, alphahat, R, RA, N, Nref, memoize=False):
        Rajd = Nadjust_after = None
        if self.params.Radjust == 'after':
            Nadjust_after = Nref
            Radj = R
        elif self.params.Radjust == 'before':
            Nadjust_after = None
            Radj = R.adjusted_before_inversion(Nref)
        else:
            Nadjust_after = None
            Radj = R

        if self.params.RAreg:
            print('regularizing RA')
            RA = RA.add_ridge(self.params.Lambda, renormalize=True)
            gs = GenomicSubset(self.params.region)
            A = SnpSubset(self.refpanel, bedtool=gs.bedtool)
            RA.zero_outside_irs(A.irs)

        if not memoize or not hasattr(self, 'bias'):
            print('adding lambda')
            Radjreg = Radj.add_ridge(self.params.Lambda, renormalize=True)
            print('computing inverse')
            self.Radjreginv = Radjreg.inv(Nadjust_after=Nadjust_after)

            print('done.computing bias...')
            A = SnpSubset(self.refpanel, bedtool=GenomicSubset(self.params.region).bedtool)
            W = self.window(A)
            if not self.params.avgunbiased:
                tr = self.Radjreginv.dot(RA).trace()
                self.scaling = 1
            else:
                tr = RA.dot(self.Radjreginv).dot(R).dot(self.Radjreginv).trace()
                Q = R.dot(self.Radjreginv).dot(RA).dot(self.Radjreginv).dot(R)
                Q.zero_outside_irs(A.irs)
                self.scaling = A.num_snps() / Q.trace()
            # self.bias = tr / N + \
            #         float(self.refpanel.M-len(W.irs))/self.refpanel.M * \
            #             self.params.sigma2g * tr / self.params.pop_size
            self.bias = tr / N + \
                        self.params.sigma2g * tr / self.params.pop_size
            print('\nbias =', self.bias)
            print('scaling =', self.scaling)

        betahat = self.Radjreginv.dot(alphahat)

        return self.scaling * (betahat.dot(RA.dot(betahat)) - self.bias)


if __name__ == '__main__':
    import primitives.genome, primitives.simulation, primitives.dataset
    reload(primitives.genome)
    reload(primitives.simulation)
    reload(primitives.dataset)
    from primitives.genome import GenomicSubset, SnpSubset
    from primitives.simulation import SumstatSimulation
    from primitives.dataset import Dataset
    est = MLE_reg_noband(refpanel='GERA_ref14k', region='50', ld_window=1, units='cM', Lambda=0.05)
    # est = MLE_reg_noband(refpanel='tinyGERA_ref2k', region='tiny', ld_window=0.01,
            # Lambda=0.01)
    sim = SumstatSimulation('GERA.50bigger_sparse')
    # est.preprocess()
    est.run(1, sim)
