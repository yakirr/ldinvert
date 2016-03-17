from __future__ import print_function, division
import numpy as np
import pandas as pd
from pybedtools import BedTool
from primitives import Dataset, GenomicSubset, SnpSubset
import paths
import argparse
import time


def compute_statistic(X, alphahat, N, A, Nref, numpcs=0):
    t0 = time.time()
    print(time.time()-t0, 'computing R')
    R = X.T.dot(X) / Nref
    print(np.count_nonzero(R), 'nonzero entries in R')
    RA = np.copy(R)
    RA[np.logical_not(A)] = 0
    RA.T[np.logical_not(A)] = 0

    print(time.time()-t0, 'computing Rri')
    Rri = np.linalg.inv(((Nref/(Nref-len(R)-1))*R + 0.03 * np.eye(len(R)))/1.03)
    # U, svs, VT = np.linalg.svd(R)
    # Rri = U[:,:numpcs].dot(np.diag(1/svs[:numpcs])).dot(VT[:numpcs,:])
    Z = Rri.dot(RA).dot(Rri)
    Q = R.dot(Z).dot(R)
    scaling = np.sum(A) / np.trace(Q)
    # scaling = 1
    print(time.time()-t0, 'scaling', scaling)
    betahat = Rri.dot(alphahat)
    point_estimate = scaling * (betahat.dot(RA).dot(betahat) - np.sum(RA * Rri) / N)
    print(time.time()-t0, 'estimate', point_estimate)

    def term1_coeff():
        return 4*(1/N) + 4/N*(1/N)
    def term2_coeff():
        return 4/N + 2/N**2
    def term3_coeff():
        return 2*(1/N)**2

    QZ = Q.dot(Z)
    QZR = QZ.dot(R)

    # term A: beta^T RZRZR beta = beta^T QZR beta
    variance1 = max(0,scaling**2 * point_estimate / np.sum(A) * np.trace(QZR) * term1_coeff())
    # term B: (beta^T Q beta)^2
    variance2 = scaling**2 * (point_estimate / scaling)**2 * term2_coeff()
    # compute the term that doesn't depend on beta
    variance3 = scaling**2 * np.trace(QZ) * term3_coeff()
    variance = variance1 + variance2 + variance3
    print(time.time()-t0, 'variance is {} + {} + {} = {}'.format(variance1, variance2, variance3, variance))
    print('zscore:', point_estimate / np.sqrt(variance))
    return point_estimate, variance, R, RA

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--refpanel', type=str, required=True)
    parser.add_argument('--ldblocks', type=str, required=False,
            default='pickrell_ldblocks.hg19.eur.bed')
    parser.add_argument('--region', type=str, required=True)
    parser.add_argument('--sumstats_path', type=str, required=True)

    args = parser.parse_args()

    print('loading reference panel')
    refpanel = Dataset(args.refpanel)

    print('loading region')
    A = GenomicSubset(args.region)

    print('loading ld blocks')
    blocks = BedTool(paths.reference + args.ldblocks)

    print('finding ld blocks that overlap with A')
    relevant_blocks = blocks.intersect(A.bedtool, wa=True).saveas()
    print('found', len(relevant_blocks), 'blocks that overlap with A')

    print('reading refpanel bim')
    refpanel_bim = pd.read_csv(refpanel.genotypes_bedfile.filename + '.bim',
            sep='\t',
            names=['CHR', 'SNP', 'cM', 'BP', 'A1', 'A2'])
    refpanel_bim['INDEX'] = np.arange(len(refpanel_bim))
    refpanel_bim['A'] = 1
    refpanel_bim.ix[SnpSubset(refpanel, A.bedtool).irs,'A'] = 1

    print('reading sumstats')
    sumstats = pd.read_csv(args.sumstats_path+'.gz', header=0, sep='\t',
            compression='gzip')

    for block in relevant_blocks:
        window_ref = SnpSubset(refpanel, BedTool([block]))
        refpanel_bim_w = refpanel_bim.loc[window_ref.irs]

        print('merging')
        refpanel_with_sumstats = refpanel_bim_w.merge(sumstats, how='left', on=['SNP'])
        refpanel_to_use = refpanel_with_sumstats.loc[
                refpanel_with_sumstats['N'].notnull()]
        refpanel_to_use = refpanel_to_use.sort(['BP'])

        print('flipping snps to line up with refpanel')
        flip = refpanel_to_use['A1_x'] == refpanel_to_use['A2_y']
        refpanel_to_use.ix[flip, 'Z'] *= -1

        print('reading in data')
        X = refpanel.get_standardized_genotypes_in_iter(refpanel_to_use['INDEX'])
        alphahat = (refpanel_to_use['Z'] / np.sqrt(refpanel_to_use['N'])).values
        N = np.mean(refpanel_to_use['N'])
        A = refpanel_to_use['A'].values

        # snps = np.random.choice(X.shape[1], size=645, replace=False)
        snps = np.sort(np.random.choice(X.shape[1], size=min(5*645,X.shape[1]), replace=False))
        X = X[:,snps]
        alphahat = alphahat[snps]
        A = A[snps]
        break

        estimate, variance, R, RA = compute_statistic(
                X, alphahat, N, A, X.shape[0])
