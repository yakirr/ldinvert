from __future__ import print_function, division
import argparse
import pandas as pd
import numpy as np
import gzip
from primitives import Dataset, SnpSubset, GenomicSubset
from pyutils import bsub
import paths
import pickle


def results_filename(annot_stem, sumstats_stem, second_annot_stem=None, chrnum=None):
    annot_name = annot_stem.split('/')[-1]
    if second_annot_stem is not None:
        second_annot_name = '-' + second_annot_stem.split('/')[-1]
    else:
        second_annot_name = ''
    pheno_name = sumstats_stem.split('/')[-1]
    stem = '/groups/price/yakir/results/2016.03.21_signed/' + \
            pheno_name + '.' + annot_name + second_annot_name
    if chrnum is None:
        return stem
    else:
        return stem + '.' + str(chrnum)

def main(args):
    print('chrom', args.chrnum)
    print('reading sumstats and filtering by sample size')
    sumstats = pd.read_csv(args.sumstats_stem + '.sumstats.gz', header=0, sep='\t',
            compression='gzip')
    N90 = np.percentile(sumstats['N'], 90)
    Nthresh = 0.9 * N90
    print('\tthreshold sample size:', Nthresh)
    print('\toriginally at', len(sumstats), 'SNPs')
    sumstats = sumstats.loc[
            sumstats['N'] >= Nthresh]
    print('\tafter filtering by N, now at', len(sumstats), 'SNPs')

    print('reading annot file')
    annot = pd.read_csv('{}.{}.annot.gz'.format(args.annot_stem, args.chrnum),
            header=0, sep='\t', compression='gzip')
    if args.second_annot_stem is not None:
        annot2 = pd.read_csv('{}.{}.annot.gz'.format(args.second_annot_stem, args.chrnum),
                header=0, sep='\t', compression='gzip')
        annot2_names = annot2.columns[4:].values.tolist()
        annot2 = annot2[['SNP'] + annot2_names]
        annot = annot.merge(annot2, how='left', on='SNP')
    annot_names = annot.columns[4:]
    print('\tannotation names:', annot_names.values)
    print('\tannotations supported on',
            np.sum(annot[annot_names].values != 0, axis=0), 'SNPs')

    print('reading refpanel bim')
    refpanel = Dataset(args.refpanel + '.' + str(args.chrnum))
    refpanel_bim = pd.read_csv(refpanel.genotypes_bedfile.filename + '.bim',
            sep='\t',
            names=['CHR', 'SNP', 'cM', 'BP', 'A1', 'A2'])
    refpanel_bim['INDEX'] = np.arange(len(refpanel_bim))
    print('\trefpanel contains', len(refpanel_bim), 'SNPs')

    print('merging files')
    data = refpanel_bim.merge(sumstats, how='left', on=['SNP'])
    data = data.loc[
            data['N'].notnull()]
    data = data.merge(annot, how='left', on=['SNP'])
    print('\tnumber of SNPs in common between refpanel and sumstats is', len(data))
    print('\tannotations now supported on',
            np.sum(data[annot_names].values != 0, axis=0), 'SNPs')

    print('\tflipping sumstats to line up with refpanel')
    flip = data['A1_x'] == data['A2_y']
    data.ix[flip, 'Z'] *= -1

    N = np.min(data['N'].values)
    alphahat = data['Z'].values / np.sqrt(data['N'].values)
    V = data[annot_names].values
    indices = data['INDEX'].values
    print('N=', N, 'alphahat.shape=', alphahat.shape, 'V.shape=', V.shape,
            'indices.shape=', indices.shape)

    def chunker(seq, size):
        return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))
    # print('computing XO')
    # XO = np.zeros(refpanel.N)
    # for group in chunker(range(len(data)), 10000):
    #     print('\tXO', group[0], group[-1])
    #     X = refpanel.get_standardized_genotypes_in_iter(indices[group])
    #     XO += X.dot(np.ones(len(group)))

    # print('computing XTXO')
    # XTXO = np.zeros(V.shape[0])
    # for group in chunker(range(len(data)), 10000):
    #     print('\tXTXO', group[0], group[-1])
    #     X = refpanel.get_standardized_genotypes_in_iter(indices[group])
    #     XTXO[group] = X.T.dot(XO)
    # RO = XTXO / refpanel.N
    # wV = V / np.abs(RO)[:,None]
    wV = V

    print('computing XV')
    XV = np.zeros((refpanel.N, V.shape[1]))
    XwV = np.zeros((refpanel.N, wV.shape[1]))
    for group in chunker(range(len(data)), 10000):
        print('\tXV', group[0], group[-1])
        X = refpanel.get_standardized_genotypes_in_iter(indices[group])
        XV += X.dot(V[group])
        XwV += X.dot(wV[group])

    print('computing XTXV')
    XTXV = np.zeros(V.shape)
    XTXwV = np.zeros(wV.shape)
    for group in chunker(range(len(data)), 10000):
        print('\tXTXV', group[0], group[-1])
        X = refpanel.get_standardized_genotypes_in_iter(indices[group])
        XTXV[group] = X.T.dot(XV)
        XTXwV[group] = X.T.dot(XwV)
    RV = XTXV / refpanel.N
    RwV = XTXwV / refpanel.N

    print('computing wV^TRV')
    wVTRV = wV.T.dot(RV)
    wVTRwV = wV.T.dot(RwV)
    print('\twV^TRV =', wVTRV)

    wVTalphahat = wV.T.dot(alphahat)
    print('chr', args.chrnum, ':', np.linalg.solve(wVTRV, wVTalphahat))

    pickle.dump((wVTalphahat, wVTRV, wVTRwV, N),
            open(results_filename(
                args.annot_stem, args.sumstats_stem,
                second_annot_stem=args.second_annot_stem,
                chrnum=args.chrnum), 'w'))

def submit(args):
    my_args = ['--annot_stem', args.annot_stem,
        '--sumstats_stem', args.sumstats_stem,
        '--refpanel', args.refpanel] + \
        ([] if args.second_annot_stem is None else ['--second_annot',args.second_annot_stem])+\
        ['main',
        '--chrnum', '$LSB_JOBINDEX']
    outfilepath = \
        results_filename(args.annot_stem, args.sumstats_stem,
                second_annot_stem=args.second_annot_stem,
                chrnum='%I') + '.log'

    bsub.submit(
            ['python', '-u', __file__] + my_args,
            outfilepath,
            jobname='run[1-22]', debug=args.debug)

def merge(args):
    wVTalphahat = 0
    wVTRV = 0
    wVTRwV = 0
    # ssq = 0
    N = 0
    for chrnum in range(1,23):
        if chrnum == 6:
            continue
        my_wVTalphahat, my_wVTRV, my_wVTRwV, N = pickle.load(open(
            results_filename(
                args.annot_stem, args.sumstats_stem,
                second_annot_stem=args.second_annot_stem,
                chrnum=chrnum)))
        print('chr {}: est={}, wVTRV={}'.format(
            chrnum,
            np.linalg.solve(my_wVTRV, my_wVTalphahat),
            my_wVTRV))
        wVTalphahat += my_wVTalphahat
        wVTRV += my_wVTRV
        wVTRwV += my_wVTRwV
        # ssq += my_wVTalphahat.T.dot(my_wVTalphahat)
    # corr = wVTalphahat / np.sqrt(wVTRV * args.h2g)
    mu = np.linalg.solve(wVTRV, wVTalphahat)
    print('wVTRV=', wVTRV)
    wVTRVi = np.linalg.inv(wVTRV)
    # variance = 1/N*wVTRVi + 1/N * wVTRVi.dot(ssq).dot(wVTRVi)
    variance = 1/N*wVTRVi.dot(wVTRwV).dot(wVTRVi.T)
    print('covariance of bethat=', variance)

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
    with open(results_filename(args.annot_stem, args.sumstats_stem,
        second_annot_stem=args.second_annot_stem), 'w') as f:
        f.write(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--annot_stem', type=str, required=True,
            help='path to annot file, not including chromosome number and extension.  ' + \
                    'For example: path/to/annot')
    parser.add_argument('--second_annot_stem', type=str, required=False, default=None)
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
    merge_parser.add_argument('--h2g', type=float, required=True,
            help='the value of h2g to use for computing genetic correlation')

    bsub.choose_parser_and_run(parser)

