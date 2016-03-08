from __future__ import print_function, division
import numpy as np
import math
import argparse
import cPickle as pickle
import os, time
from pybedtools import BedTool
from pysnptools.util import IntRangeSet
from estimator import Estimator
from primitives import Dataset, GenomicSubset, SnpSubset, SnpPartition
from sparse.blockdiag import BlockDiag
import paths
from pyutils import fs

class Ldinvert(Estimator):
    parser = argparse.ArgumentParser(add_help=False, parents=[Estimator.parser])
    parser.add_argument('--breakpointsfile', type=str, required=True,
            help='the name of the LD breakpoints bed file to use. This file should be \
                    in data/reference')
    parser.add_argument('--region', type=str, required=True,
            help='the name of the subset of the genome whose heritability should be \
                    analyzed. these files are in data/genome_subsets')
    parser.add_argument('--Lambda', type=float, required=True,
            help='the value of lambda to use in regularization')
    parser.add_argument('--pop_size', type=int, required=False, default=10**9,
            help='the population size to adjust for when accounting for effects of snps \
                    outside the window. Defaults to essentially no correction.')
    parser.add_argument('--num_chunks', type=int, required=False, default=1,
            help='the number of chunks the LD blocks are getting split into')
    parser.add_argument('--chunk', type=int, required=False, default=1,
            help='the 1-based index of the chunk to analyze')
    parser.add_argument('--preprocessonfly', type=int, required=False, default=0,
            help='0 means do not preprocess on the fly; 1 means do the preprocessing' + \
                    'on the fly')

    def readable_name(self):
        return 'Ldinvert,A={},L={:0.3f},ref={},bpf={},Npop={}'.format(
                self.params.region,
                self.params.Lambda,
                self.params.refpanel,
                self.params.breakpointsfile,
                self.params.pop_size)

    def covariance_path(self):
        return self.refpanel.auxfiles_path + 'pre.cov.bpf={}/'.format(
                self.params.breakpointsfile)
    def covariance_preprocessing_in_progress(self):
        return os.path.exists(self.covariance_path())
    def declare_covariance_preprocessing_in_progress(self):
        fs.makedir(self.covariance_path())

    def invcovariance_path(self):
        return self.refpanel.auxfiles_path + 'pre.invcov.L={:0.3f}.bpf={}/'.format(
                self.params.Lambda,
                self.params.breakpointsfile)
    def invcovariance_preprocessing_in_progress(self):
        return os.path.exists(self.invcovariance_path())
    def declare_invcovariance_preprocessing_in_progress(self):
        fs.makedir(self.invcovariance_path())

    def preprocessing_foldername(self):
        return 'pre.ldinvert.A={}.L={:0.3f}.bpf={}'.format(
                self.params.region,
                self.params.Lambda,
                self.params.breakpointsfile)

    def R_file(self, mode='rb'):
        return open(self.covariance_path() + 'R.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks), mode)
    def Rri_file(self, mode='rb'):
        return open(self.invcovariance_path() + 'Rri.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks), mode)
    def RA_file(self, mode='rb'):
        return open(self.path_to_preprocessed_data() + 'RA.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks), mode)
    def biasmatrix_file(self, mode='rb'):
        return open(self.path_to_preprocessed_data() + 'ZR.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks), mode)
    def set_scalings(self, scalings):
        pickle.dump(scalings, open(self.path_to_preprocessed_data() + 'scaling.{}of{}'.format(
            self.params.chunk,
            self.params.num_chunks), 'w'))
    def get_scalings(self):
        return pickle.load(open(self.path_to_preprocessed_data() + 'scaling.{}of{}'.format(
            self.params.chunk,
            self.params.num_chunks), 'r'))
    def save_variance_matrices(self, Q, Z, QZ, QZR):
        pickle.dump(Q, open(self.path_to_preprocessed_data() + 'Q.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks), 'wb'), 2)
        pickle.dump(Z, open(self.path_to_preprocessed_data() + 'Z.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks), 'wb'), 2)
        pickle.dump(QZ, open(self.path_to_preprocessed_data() + 'QZ.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks), 'wb'), 2)
        pickle.dump(QZR, open(self.path_to_preprocessed_data() + 'QZR.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks), 'wb'), 2)
    def get_variance_matrices(self):
        return (pickle.load(open(self.path_to_preprocessed_data() + 'Q.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks))),
            pickle.load(open(self.path_to_preprocessed_data() + 'Z.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks))),
            pickle.load(open(self.path_to_preprocessed_data() + 'QZ.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks))),
            pickle.load(open(self.path_to_preprocessed_data() + 'QZR.{}of{}.bd'.format(
            self.params.chunk,
            self.params.num_chunks))))
    def init(self):
        self.Rri = pickle.load(self.Rri_file())
        self.R = pickle.load(self.R_file())
        self.RA = pickle.load(self.RA_file())
        self.A = SnpSubset(self.refpanel, GenomicSubset(self.params.region).bedtool)
        self.ZR = pickle.load(self.biasmatrix_file())
        self.Q, self.Z, self.QZ, self.QZR = self.get_variance_matrices()

    def chunk_size(self, ranges):
        return int(math.ceil(len(ranges) / self.params.num_chunks))
    def ranges_in_chunk(self, ranges):
        start = (self.params.chunk-1)*self.chunk_size(ranges)
        end = min(len(ranges), start + self.chunk_size(ranges))
        return ranges[start:end]

    def preprocess_memoryreq_GB(self):
        return 8
    def preprocess(self, use_filesystem=True):
        if not self.covariance_preprocessing_in_progress():
            print('covariance matrix not found. creating...')
            self.declare_covariance_preprocessing_in_progress()
            self.compute_covariance()
        print('loading covariance matrix')
        R = pickle.load(self.R_file())

        if not self.invcovariance_preprocessing_in_progress():
            print('inverse covariance matrix not found. creating...')
            self.declare_invcovariance_preprocessing_in_progress()
            self.compute_invcovariance()
        print('loading inverse covariance matrix')
        Rri = pickle.load(self.Rri_file())

        t0 = time.time()
        print(time.time() - t0, ': creating and saving RA')
        A = SnpSubset(self.refpanel, GenomicSubset(self.params.region).bedtool)
        RA = R.copy(); RA.zero_outside_irs(A.irs)
        pickle.dump(RA, self.RA_file(mode='wb'), 2)

        print(time.time() - t0, ': computing and saving scaling')
        Z = Rri.dot(RA.dot(Rri))
        Q = R.dot(Z).dot(R)
        QA = Q.copy()
        QA.zero_outside_irs(A.irs)
        scalings = {r :
                len(A.irs & IntRangeSet(r)) / np.trace(QA.ranges_to_arrays[r])
                for r in QA.ranges()}
        print(time.time() - t0, ': scalings are', scalings)
        self.set_scalings(scalings)

        print(time.time() - t0, ': computing and saving bias matrix')
        ZR = RA.dot(Rri).dot(R).dot(Rri)
        pickle.dump(ZR, self.biasmatrix_file(mode='wb'), 2)

        print(time.time() - t0, ': variance matrices')
        QZ = Q.dot(Z)
        QZR = QZ.dot(R)
        self.save_variance_matrices(Q, Z, QZ, QZR)
        print(time.time() - t0, ': done')

    def compute_covariance(self):
        breakpoints = BedTool(paths.reference + self.params.breakpointsfile)
        blocks = SnpPartition(self.refpanel, breakpoints, remove_mhc=True)
        R = BlockDiag.ld_matrix_blocks(self.refpanel, self.ranges_in_chunk(blocks.ranges()))
        pickle.dump(R, self.R_file(mode='wb'), 2)

    def compute_invcovariance(self):
        R = pickle.load(self.R_file())
        if '_ref' in self.refpanel.name:
            # print('this is a reference panel, so Im adjusting for the fact that the sample'+\
                    # ' precision matrix will be inflated')
            # R = R.adjusted_before_inversion(self.refpanel.N)
            pass
        Rreg = R.add_ridge(self.params.Lambda, renormalize=True)
        Rri = Rreg.inv(Nadjust_after=None)
        pickle.dump(Rri, self.Rri_file(mode='wb'), 2)

    # definitions:
    # Z = Rri RA Rri
    # Q = RZR
    # expectation of alphahat Rri RA Rri alphahat = alphahat Z alphahat is
    #   beta Q beta + Tr(ZR)/N
    def run(self, beta_num, sim):
        if not hasattr(self, 'R'):
            print('loading matrices')
            self.init()
        self.beta = pickle.load(sim.beta_file(beta_num))

        print('computing bias')
        self.biases = self.compute_biases(sim.sample_size)
        print('biases are', self.biases)
        self.scalings = self.get_scalings()

        # compute the results
        results = []
        variances = []
        for alphahat in sim.sumstats_aligned_to_refpanel(beta_num, self.refpanel):
            alphahat = BlockDiag.from_big1darray(alphahat, self.R.ranges())
            results.append(self.compute_statistic(
                alphahat))
            variances.append(self.compute_variance(
                alphahat, results[-1], sim.sample_size))
            print(len(results), results[-1], variances[-1])

        print('empirical var of results:', np.var(results))

        return np.concatenate([np.array([results]).T,
            np.array([variances]).T], axis=1)

    def c(self, r):
        return np.linalg.norm(np.concatenate(
            [self.beta[:r[0]], self.beta[r[1]:]]))**2 / self.params.pop_size

    def compute_biases(self, N):
        M = self.refpanel.M
        biases = { r:
            np.trace(W) * (1/N + self.c(r))
            for r, W in self.ZR.ranges_to_arrays.items()}
        return biases

    def compute_statistic(self, alphahat):
        betahat = self.Rri.dot(alphahat)
        return sum([
            self.scalings[r] * (betahat.ranges_to_arrays[r].dot(
                self.RA.ranges_to_arrays[r]).dot(
                betahat.ranges_to_arrays[r]) - \
                        self.biases[r])
            for r in self.RA.ranges()])

    def compute_variance(self, alphahat, point_estimate, N, use_beta=False):
        def term1_coeff(r):
            return 4*(1/N+self.c(r)) + 4/N*(1/N+self.c(r))
        def term2_coeff():
            return 4/N + 2/N**2
        def term3_coeff(r):
            return 2*(1/N+self.c(r))**2

        # compute the term that doesn't depend on beta
        variance3 = sum([
            self.scalings[r]**2 * np.trace(self.QZ.ranges_to_arrays[r]) * \
                    term3_coeff(r)
            for r in self.Q.ranges()])

        # now compute the other two terms
        if use_beta:
            beta = self.beta
            # term A: beta^T RZRZR beta = beta^T QZR beta
            variance1 = sum([
                self.scalings[r]**2 * \
                    beta[r[0]:r[1]].dot(self.QZR.ranges_to_arrays[r].dot(beta[r[0]:r[1]])) * \
                    term1_coeff(r)
                for r in self.Q.ranges()])
            # term B: (beta^T Q beta)^2
            variance2 = sum([
                self.scalings[r]**2 * \
                    beta[r[0]:r[1]].dot(self.Q.ranges_to_arrays[r].dot(beta[r[0]:r[1]]))**2 * \
                    term2_coeff()
                for r in self.Q.ranges()])
        else:
            # term A
            # alphahatTZR = alphahat.dot(self.ZR)
            # Zalphahat = self.Z.dot(alphahat)
            # termAbiases = {r :
            #     np.einsum('ij,ji',self.ZR.ranges_to_arrays[r],self.ZR.ranges_to_arrays[r]) * \
            #             (1/N + self.c(r))
            #     for r in self.Q.ranges()}
            # variance1 = sum([
            #     (alphahatTZR.ranges_to_arrays[r].dot(Zalphahat.ranges_to_arrays[r]) - \
            #             termAbiases[r]) * \
            #             term1_coeff(r)
            #     for r in self.Q.ranges()])
            betahat = self.Rri.dot(alphahat)
            point_estimates = {r:
                    self.scalings[r] * (betahat.ranges_to_arrays[r].dot(
                        self.RA.ranges_to_arrays[r]).dot(
                        betahat.ranges_to_arrays[r]) - self.biases[r])
                for r in self.Q.ranges()}
            variance1 = sum([
                self.scalings[r]**2 * point_estimates[r] / len(self.A.irs & IntRangeSet(r)) * \
                    np.trace(self.QZR.ranges_to_arrays[r]) * \
                    term1_coeff(r)
                for r in self.Q.ranges()])

            # term B
            variance2 = term2_coeff() * sum([
                self.scalings[r]**2 * (point_estimates[r] / self.scalings[r])**2
                for r in self.Q.ranges()])

        variance = variance1 + variance2 + variance3
        print('\nvariance is {} + {} + {} = {}'.format(variance1, variance2, variance3, variance))
        return variance
