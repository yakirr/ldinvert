from __future__ import print_function, division
import argparse
import pickle
import pandas as pd
import numpy as np
import gzip
from pybedtools import BedTool
from primitives import Dataset, SnpSubset, GenomicSubset, SnpPartition
from pyutils import bsub
import paths
import signed_fe


def results_filename(annot_stems, sumstats_stem, chrnum=None):
    return signed_fe.results_filename(annot_stems, sumstats_stem, chrnum=chrnum) + '.re'

def mult_by_R(V, data, refpanel, indices):
    def chunker(seq, size):
        return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))
    print('computing XV')
    XV = np.zeros((refpanel.N, V.shape[1]))
    for group in chunker(range(refpanel.M), 10000):
        print('\tXV', group[0], group[-1])
        X = refpanel.get_standardized_genotypes_in_iter(indices[group])
        XV += X.dot(V[group])

    print('computing XTXV')
    XTXV = np.zeros(V.shape)
    for group in chunker(range(refpanel.M), 10000):
        print('\tXTXV', group[0], group[-1])
        X = refpanel.get_standardized_genotypes_in_iter(indices[group])
        XTXV[group] = X.T.dot(XV)

    return XTXV / refpanel.N

def get_ldscores(args):
    print('reading ld scores')
    ldscfilename = '/groups/price/ldsc/reference_files/' + \
            'UK10K/totalld_allsnps/ALSPAC.TWINSUK.QC.' + str(args.chrnum) + '.l2.ldscore.gz'
    ldscores = pd.read_csv(ldscfilename, header=0, sep='\t', compression='gzip')
    return ldscores[['SNP', 'L2']]

def main(args):
    print('chrom', args.chrnum)
    sumstats = signed_fe.get_sumstats(args)
    annot, annot_names = signed_fe.get_annot(args)
    refpanel, refpanel_bim = signed_fe.get_refpanel(args)
    ldscores = get_ldscores(args)
    N, alphahat, V, ss_indices, l2_ss = signed_fe.merge_files(
            refpanel_bim, sumstats, annot, annot_names, get_full_annot=True, ldscores=ldscores)

    # compute weights
    print('computing weights, assuming h2g=0.25')
    sigma2g = 0.25 / refpanel.M
    weights = l2_ss / N + sigma2g * l2_ss**2
    print('\tmean weight: {}, stddev: {}'.format(np.mean(1./weights), np.std(1./weights)))

    RV = signed_fe.mult_by_R_ldblocks(V, refpanel)

    print('computing estimate and covariance')
    # get only regression snps
    RV_ss = RV[ss_indices]
    # weight each sample by the right weight
    RV_ss_w = RV_ss / weights[:,None]
    PiT_RV_ss_w = np.zeros((refpanel.M, RV_ss_w.shape[1]))
    PiT_RV_ss_w[ss_indices] = RV_ss_w
    inter1 = signed_fe.mult_by_R_ldblocks(PiT_RV_ss_w, refpanel)
    inter2 = signed_fe.mult_by_R_ldblocks(inter1, refpanel)
    M1 = PiT_RV_ss_w.T.dot(inter1)
    M2 = PiT_RV_ss_w.T.dot(inter2)

    VTRalphahat = RV_ss_w.T.dot(alphahat)
    VTR2V =  RV_ss_w.T.dot(RV_ss)
    print('VTR2V=', VTR2V)
    point_estimate = np.linalg.solve(VTR2V, VTRalphahat)
    print('chr', args.chrnum, ':', point_estimate)

    pickle.dump((VTRalphahat, VTR2V, M1, M2, N),
            open(results_filename(
                args.annot_stems, args.sumstats_stem,
                chrnum=args.chrnum), 'w'))

def submit(args):
    my_args = ['--annot_stems', ' '.join(args.annot_stems),
        '--sumstats_stem', args.sumstats_stem,
        '--refpanel', args.refpanel,
        'main',
        '--chrnum', '$LSB_JOBINDEX']
    outfilepath = \
        results_filename(args.annot_stems, args.sumstats_stem,
                chrnum='%I') + '.log'

    bsub.submit(
            ['python', '-u', __file__] + my_args,
            outfilepath,
            jobname='run[1-22]', debug=args.debug)

def merge(args):
    VTRalphahat = 0
    VTR2V = 0
    M1 = 0
    M2 = 0
    for chrnum in range(1,23):
        my_VTRalphahat, my_VTR2V, my_M1, my_M2, N = pickle.load(open(
            results_filename(args.annot_stems, args.sumstats_stem,
                chrnum=chrnum)))
        print('chr {}: {}, VTR2V={}'.format(
            chrnum,
            my_VTRalphahat,
            my_VTR2V))
        VTRalphahat += my_VTRalphahat
        VTR2V += my_VTR2V
        M1 += my_M1
        M2 += my_M2
    mu = np.linalg.solve(VTR2V, VTRalphahat)
    print('mu=', mu)

    variance1 = np.linalg.solve(VTR2V, M1)
    variance1 = np.linalg.solve(VTR2V.T, variance1.T).T / N
    variance2 = np.linalg.solve(VTR2V, M2)
    variance2 = np.linalg.solve(VTR2V.T, variance2.T).T * 0.25 / 13000000
    variance = variance1 + variance2
    print('covariance of muhat=', variance)

    output = ''
    for i, c in enumerate(mu):
        output += '{} ({}), Z={}\n'.format(
            c,
            np.sqrt(variance[i,i]),
            c / np.sqrt(variance[i,i]))
    output += 'sum: {} ({}), Z={}\n'.format(
            np.sum(mu),
            np.sqrt(np.sum(variance)),
            np.sum(mu) / np.sqrt(np.sum(variance)))
    print(output)
    with open(results_filename(args.annot_stems, args.sumstats_stem), 'w') as f:
        f.write(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--annot_stems', type=str, required=True, nargs='+',
            help='path to annot files, not including chromosome number and extension.  ' + \
                    'For example: path/to/annot path/to/annot2')
    parser.add_argument('--sumstats_stem', type=str, required=True,
            help='path to sumstats file, not including .sumstats.gz extension.')
    parser.add_argument('--refpanel', type=str, required=True,
            help='the name of the reference panel, for synchronizing allele coding')

    main_parser, submit_parser, merge_parser = bsub.add_main_and_submit(parser,
            main,
            submit,
            merge)
    main_parser.add_argument('--chrnum', type=int, required=True,
            help='which chromosome to analyze')
    submit_parser.add_argument('-debug', action='store_true', default=False)

    bsub.choose_parser_and_run(parser)

