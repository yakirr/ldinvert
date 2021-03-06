from __future__ import print_function, division
import argparse
import pandas as pd
import numpy as np
import gzip
from pybedtools import BedTool
from primitives import Dataset, SnpSubset, GenomicSubset, SnpPartition
from pyutils import bsub
import paths
import pickle


def results_filename(annot_stems, sumstats_stem, chrnum=None):
    annot_name = '-'.join([annot_stem.split('/')[-1] for annot_stem in annot_stems])
    pheno_name = sumstats_stem.split('/')[-1]
    stem = '/groups/price/yakir/results/2016.03.29_signed/' + \
            pheno_name + '.' + annot_name
    if chrnum is None:
        return stem
    else:
        return stem + '.' + str(chrnum)

def get_sumstats(args):
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
    return sumstats

def get_annot(args):
    print('reading annot files')
    annot = None
    for annot_stem in args.annot_stems:
        print('\t', annot_stem)
        if annot is None:
            annot = pd.read_csv('{}.{}.annot.gz'.format(annot_stem, args.chrnum),
                    header=0, sep='\t', compression='gzip')
        else:
            newannot = pd.read_csv('{}.{}.annot.gz'.format(annot_stem, args.chrnum),
                    header=0, sep='\t', compression='gzip')
            newannot_names = newannot.columns[4:].values.tolist()
            newannot = newannot[['SNP'] + newannot_names]
            annot = annot.merge(newannot, how='left', on='SNP')
    annot_names = annot.columns[4:]
    print('\tannotation names:', annot_names.values)
    print('\tannotations supported on',
            np.sum(annot[annot_names].values != 0, axis=0), 'SNPs')
    return annot, annot_names

def get_refpanel(args):
    print('reading refpanel and refpanel bim')
    refpanel = Dataset(args.refpanel + '.' + str(args.chrnum))
    refpanel_bim = pd.read_csv(refpanel.genotypes_bedfile.filename + '.bim',
            sep='\t',
            names=['CHR', 'SNP', 'cM', 'BP', 'A1', 'A2'])
    refpanel_bim['refpanelINDEX'] = np.arange(len(refpanel_bim))
    print('\trefpanel contains', len(refpanel_bim), 'SNPs')
    return refpanel, refpanel_bim

def merge_files(refpanel_bim, sumstats, annot, annot_names, ldscores=None,
        get_full_annot=False):
    print('merging files')
    print('\tmerging refpanel and sumstats')
    data = refpanel_bim.merge(sumstats, how='left', on=['SNP'])
    print('\tmerging annot')
    data = data.merge(annot, how='left', on=['SNP'])
    if get_full_annot:
        V = data[annot_names].values

    data = data.loc[
            data['N'].notnull()]
    print('\tnumber of SNPs in common between refpanel and sumstats is', len(data))
    print('\tannotations now supported on',
            np.sum(data[annot_names].values != 0, axis=0), 'SNPs')

    if ldscores is not None:
        print('\tmerging ldscores')
        data = data.merge(ldscores, how='left', on=['SNP'])

    print('\tflipping sumstats to line up with refpanel')
    flip = data['A1_x'] == data['A2_y']
    data.ix[flip, 'Z'] *= -1

    print('\tcreating numpy arrays')
    N = np.min(data['N'].values)
    alphahat = data['Z'].values / np.sqrt(data['N'].values)
    l2_ss = None if ldscores is None else data['L2'].values
    if not get_full_annot:
        V = data[annot_names].values
    refpanel_indices = data['refpanelINDEX'].values
    print('\tN=', N, 'alphahat.shape=', alphahat.shape, 'V.shape=', V.shape,
            'refpanel_indices.shape=', refpanel_indices.shape)

    return N, alphahat, V, refpanel_indices, l2_ss

def mult_by_R_ldblocks(V, refpanel, refpanel_indices=None):
    breakpoints = BedTool(paths.reference + 'pickrell_breakpoints.hg19.eur.bed')
    blocks = SnpPartition(refpanel, breakpoints, remove_mhc=True)

    if refpanel_indices is not None:
        Vwithzeros = np.zeros((refpanel.M, V.shape[1]))
        Vwithzeros[refpanel_indices] = V
    else:
        Vwithzeros = V
    result = np.zeros(Vwithzeros.shape)
    for r in blocks.ranges():
        print('\tXTXV', r[0], r[1], 'of', refpanel.M)
        X = refpanel.get_standardized_genotypes(r)
        result[r[0]:r[1],:] = X.T.dot(X.dot(Vwithzeros[r[0]:r[1],:]))
    if refpanel_indices is not None:
        return result[refpanel_indices] / refpanel.N
    else:
        return result / refpanel.N

def main(args):
    print('chrom', args.chrnum)
    sumstats = get_sumstats(args)
    annot, annot_names = get_annot(args)
    refpanel, refpanel_bim = get_refpanel(args)
    N, alphahat, V, refpanel_indices, _ = merge_files(
        refpanel_bim, sumstats, annot, annot_names)
    Vw = V

    print('computing RV and RVw')
    RV = mult_by_R_ldblocks(V, refpanel, refpanel_indices)
    RVw = RV #mult_by_R_ldblocks(Vw, refpanel, refpanel_indices)

    print('computing Vw^TRV')
    VwTRV = Vw.T.dot(RV)
    VwTRVw = Vw.T.dot(RVw)
    print('\tVw^TRV =', VwTRV)

    VwTalphahat = Vw.T.dot(alphahat)
    print('chr', args.chrnum, ':', np.linalg.solve(VwTRV, VwTalphahat))

    pickle.dump((VwTalphahat, VwTRV, VwTRVw, N),
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
    VwTalphahat = 0
    VwTRV = 0
    VwTRVw = 0
    N = 0
    for chrnum in range(1,23):
        if chrnum == 6:
            continue
        my_VwTalphahat, my_VwTRV, my_VwTRVw, N = pickle.load(open(
            results_filename(
                args.annot_stems, args.sumstats_stem,
                chrnum=chrnum)))
        print('chr {}: est={}, VwTRV={}'.format(
            chrnum,
            np.linalg.solve(my_VwTRV, my_VwTalphahat),
            my_VwTRV))
        VwTalphahat += my_VwTalphahat
        VwTRV += my_VwTRV
        VwTRVw += my_VwTRVw
    mu = np.linalg.solve(VwTRV, VwTalphahat)
    print('VwTRV=', VwTRV)
    wVTRVi = np.linalg.inv(VwTRV)
    # variance = 1/N*wVTRVi + 1/N * wVTRVi.dot(ssq).dot(wVTRVi)
    variance = 1/N*wVTRVi.dot(VwTRVw).dot(wVTRVi.T)
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
            help='paths to annot files, not including chromosome number and extension.  ' + \
                    'For example: path/to/annot1 path/to/annot2')
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

