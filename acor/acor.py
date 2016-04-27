from __future__ import print_function, division
import argparse
import numpy as np
import scipy.stats as stats
import pandas as pd
import pyutils.pretty as pretty
import primitives.dataset as prd
import primitives.annotation as pa
import common


def find_conv1_names(names):
    return [n for n in names if '.conv1' in n]
def find_conv3_names(names):
    return [n for n in names if '.conv3' in n]
def make_conv1_names(names):
    return [n+'.conv1' for n in names]
def make_conv3_names(names):
    return [n+'.conv3' for n in names]

def cor_fe(args):
    print('FIXED EFFECTS')
    print('loading information by chromosome')
    annot, conv = [], []
    refpanel = prd.Dataset(args.bfile_chr)
    annots = [pa.Annotation(fname) for fname in args.annot_chr]
    convannots = [pa.Annotation(fname) for fname in args.conv_chr]

    # load annot and conv files and merge across chromosomes
    for chrnum in args.chroms:
        print('== chr', chrnum)
        # read data
        myannot, annot_names = common.get_annots(
                [a.sannot_filename(chrnum) for a in annots])
        myconv, conv_names = common.get_convs(
                [a.conv_filename(chrnum, args.fullconv) for a in convannots],
                [a.sannot_filename(chrnum) for a in annots],
                (refpanel, chrnum),
                args.ld_breakpoints,
                args.mhc_path)
        # basic error checks
        if (myannot['SNP']!=myconv['SNP']).any():
            raise Exception(
                    'ERROR: annot and conv do not contain identical snps in the same order')
        if find_conv1_names(myconv.columns) != make_conv1_names(annot_names):
            raise Exception(
                    'ERROR: conv file must contain same columns in same order as annot file')
        # append to previous chromosomes
        annot.append(myannot)
        conv.append(myconv)
    annot = pd.concat(annot, axis=0).reset_index(drop=True)
    conv = pd.concat(conv, axis=0).reset_index(drop=True)

    # merge annot and conv into one dataframe with the right columns
    annot[conv_names] = \
            conv[conv_names]
    print('==done')

    # load sumstats and merge them with the annot, flipping alleles if necessary
    sumstats = common.get_sumstats(args.sumstats)
    print('merging sumstats and annot file')
    zannotconv = pd.merge(sumstats, annot, how='inner', on='SNP')
    print('matching sumstats alleles and annot alleles')
    toflip = zannotconv['A1_x'] == zannotconv['A2_y']
    zannotconv.ix[toflip, 'Z'] *= -1

    # create the relevant numpy arrays
    N = np.min(zannotconv['N'].values)
    V = zannotconv[annot_names].values
    RV = zannotconv[conv_names].values
    alphahat = zannotconv['Z'] / np.sqrt(zannotconv['N'].values)

    # compute the estimate
    VTRV = V.T.dot(RV)
    estimate = np.linalg.solve(VTRV, V.T.dot(alphahat))
    cov = np.linalg.inv(VTRV) / N

    # print output
    output = pd.DataFrame(columns=('NAME','MU_EST','MU_STDERR','MU_Z','MU_P', 'TOP', 'BOTTOM'))
    for i, n in enumerate(annot_names):
        std = np.sqrt(cov[i,i])
        output.loc[i] = [n, estimate[i], std,
                estimate[i]/std,
                2*stats.norm.sf(estimate[i]/std,0,1),
                ','.join(map(str, LambdaRV.T.dot(alphahat).reshape((-1,)))),
                ','.join(map(str, Sigma.reshape((-1,))))]
    print(output)

    covariance = pd.DataFrame(columns=annot_names, data=cov)
    print('\nfull covariance matrix:\n{}'.format(cov))

    output.to_csv(args.out+'.results', sep='\t', index=False)
    covariance.to_csv(args.out+'.cov', sep='\t', index=False)

#TODO: replace chained calls to pd.append throughout with one call to pd.concat
def cor_re(args):
    print('RANDOM EFFECTS')

    ldscores = common.get_ldscores_allchr(args.ldscores_chr, args.chroms)
    ldscores = ldscores[['SNP','L2']]
    sumstats = common.get_sumstats(args.sumstats)
    N = np.min(sumstats['N'].values)
    print('merging sumstats and ldscores')
    sumstats = pd.merge(sumstats, ldscores, how='inner', on='SNP')
    sigma2g = (np.linalg.norm(sumstats['Z']/np.sqrt(sumstats['N']))**2 - \
            len(sumstats)/N) / np.sum(sumstats['L2'])
    sigma2g = min(max(0,sigma2g), 1/len(sumstats))
    print('h2g estimated at:', sigma2g*len(sumstats), 'sigma2g:', sigma2g)

    print('loading information by chromosome')
    annot, conv = [], []
    refpanel = prd.Dataset(args.bfile_chr)
    annots = [pa.Annotation(fname) for fname in args.annot_chr]
    convannots = [pa.Annotation(fname) for fname in args.conv_chr]

    # load annot and conv files and merge across chromosomes
    for chrnum in args.chroms:
        print('== chr', chrnum)
        # read data
        myannot, annot_names = common.get_annots(
                [a.sannot_filename(chrnum) for a in annots])
        myconv, conv_names = common.get_convs(
                [a.conv_filename(chrnum, args.fullconv) for a in convannots],
                [a.sannot_filename(chrnum) for a in annots],
                (refpanel, chrnum),
                args.ld_breakpoints,
                args.mhc_path)
        myldscores = common.get_ldscores(args.ldscores_chr + chrnum + '.l2.ldscore.gz')
        # basic error checks
        if (myannot['SNP']!=myconv['SNP']).any() or (myannot['SNP']!=myldscores['SNP']).any():
            raise Exception(
                    'ERROR: annot, conv, and ldscores do not contain identical snps ' +\
                            'in the same order')
        if find_conv1_names(myconv.columns) != make_conv1_names(annot_names):
            raise Exception(
                    'ERROR: conv file must contain same columns in same order as annot file')
        # compute weights and reconvolve
        myconv['Lambda'] = 1./(
                np.maximum(1, myldscores['L2']) / N + \
                        sigma2g * np.maximum(1, myldscores['L2'])**2)
        if args.noweights:
            print('NOT using weights')
            myconv[conv_names+'.w'] = myconv[conv_names]
        else:
            print('using weights')
            myconv[conv_names+'.w'] = myconv['Lambda'].values[:,None] * myconv[conv_names]
        myconv, _ = common.convolve(myconv, conv_names+'.w',
                (refpanel, chrnum), args.ld_breakpoints, args.mhc_path,
                fullconv=args.fullconv)
        # append to previous chromosomes
        annot.append(myannot)
        conv.append(myconv)
    annot = pd.concat(annot, axis=0).reset_index(drop=True)
    conv = pd.concat(conv, axis=0).reset_index(drop=True)

    # merge annot and conv into one dataframe with the right columns
    names = np.concatenate([conv_names,conv_names+'.w',conv_names+'.w.conv1'])
    annot[names] = \
            conv[names]
    print('==done')

    # merge annot with sumstats, flipping alleles if necessary
    print('merging sumstats and annot file')
    zannotconv = pd.merge(sumstats, annot, how='inner', on='SNP')
    print('matching sumstats alleles and annot alleles')
    toflip = zannotconv['A1_x'] == zannotconv['A2_y']
    zannotconv.ix[toflip, 'Z'] *= -1

    # create the relevant numpy arrays
    V = zannotconv[annot_names].values
    RV = zannotconv[conv_names].values
    LambdaRV = zannotconv[[n+'.w' for n in conv_names]].values
    RLambdaRV = zannotconv[[n+'.w.conv1' for n in conv_names]].values
    alphahat = zannotconv['Z'] / np.sqrt(zannotconv['N'].values)

    # compute the estimate
    Sigma = RV.T.dot(LambdaRV)
    estimate = np.linalg.solve(Sigma, LambdaRV.T.dot(alphahat))
    Sigmai = np.linalg.inv(Sigma)
    var1 = Sigmai.dot(LambdaRV.T).dot(RLambdaRV.dot(Sigmai)) / N
    var2 = np.linalg.solve(Sigma, RLambdaRV.T)
    var2 = sigma2g * var2.dot(var2.T)
    cov = var1 + var2

    # print output
    output = pd.DataFrame(columns=('NAME','MU_EST','MU_STDERR','MU_Z','MU_P', 'TOP', 'BOTTOM'))
    for i, n in enumerate(annot_names):
        std = np.sqrt(cov[i,i])
        output.loc[i] = [n, estimate[i], std,
                estimate[i]/std,
                2*stats.norm.sf(estimate[i]/std,0,1),
                ','.join(map(str, LambdaRV.T.dot(alphahat).reshape((-1,)))),
                ','.join(map(str, Sigma.reshape((-1,))))]
    print(output)

    covariance = pd.DataFrame(columns=annot_names, data=cov)
    print('\nfull covariance matrix:\n{}'.format(cov))

    output.to_csv(args.out+'.results', sep='\t', index=False)
    covariance.to_csv(args.out+'.cov', sep='\t', index=False)

def conv(args):
    print('CONV')
    print('loading information by chromosome')
    annot, conv = pd.DataFrame(), pd.DataFrame()
    refpanel = prd.Dataset(args.bfile_chr)
    annots = [pa.Annotation(fname) for fname in args.annot_chr]
    convannots = [pa.Annotation(fname) for fname in args.conv_chr]

    # load annot and conv files. they will automatically be created if they don't exist
    for chrnum in args.chroms:
        print('== chr', chrnum)
        # read data
        myconv, conv_names = common.get_convs(
                [a.conv_filename(chrnum, args.fullconv) for a in convannots],
                [a.sannot_filename(chrnum) for a in annots],
                (refpanel, chrnum),
                args.ld_breakpoints,
                args.mhc_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ld-breakpoints',
            default='/groups/price/yakir/data/reference/pickrell_breakpoints.hg19.eur.bed',
            help='path to UCSC bed file containing one zero-length bed interval per LD' + \
                    ' breakpoint')
    parser.add_argument('--mhc-path',
            default='/groups/price/yakir/data/reference/hg19.MHC.bed',
            help='path to UCSC bed file containing one zero-length bed interval per LD' + \
                    ' breakpoint')
    parser.add_argument('--annot-chr', required=True, nargs='+',
            help='space-delimited list of paths to gzipped annot files, not including ' + \
                    'chromosome number or .annot.gz extension')
    parser.add_argument('--bfile-chr', required=True,
            help='path to plink bfile of reference panel to use, not including chrom num')
    parser.add_argument('--conv-chr', nargs='+',
            help='space-delimited list of gzipped conv filenames, not including ' + \
                    'chromosome number or .conv(full).gz extension. If not provided, we use'+\
                    ' the value of --annot-chr for this')
    parser.add_argument('-fullconv', action='store_true', default=False,
            help='use/generate .fullconv.gz files, which dont take LD blocks into account')
    subparsers = parser.add_subparsers()

    subparser_conv = subparsers.add_parser('conv')
    subparser_conv.set_defaults(_func=conv)
    subparser_conv.add_argument('--chroms', nargs='+',
            default=[str(i) for i in range(1,23)],
            help='which chromosomes to analyze')

    subparser_cor = subparsers.add_parser('cor')
    subparser_cor.add_argument('--ldscores-chr', required=True,
            help='path to a set of .l2.ldscore.gz files containin a column named L2 with ' + \
                    'ld scores')
    subparser_cor.add_argument('--sumstats', required=True,
            help='path to sumstats.gz file, including extension')
    subparser_cor.add_argument('--out',
            help='stem of output files, to which we append .results and .cov')
    #TODO: add ability to input separate ld scores for weights, i.e., ld-scores TO only the
    #   regression snps.

    subparser_cor_subs = subparser_cor.add_subparsers()
    subparser_cor_fe = subparser_cor_subs.add_parser('fe')
    subparser_cor_fe.set_defaults(_func=cor_fe)
    subparser_cor_fe.add_argument('--chroms', nargs='+',
            default=[str(i) for i in range(1,23)],
            help='which chromosomes to analyze')
    subparser_cor_re = subparser_cor_subs.add_parser('re')
    subparser_cor_re.set_defaults(_func=cor_re)
    subparser_cor_re.add_argument('-noweights', action='store_true', default=False,
            help='do not weight the regression')
    subparser_cor_re.add_argument('--chroms', nargs='+',
            default=[str(i) for i in range(1,23)],
            help='which chromosomes to analyze')

    args, _ = parser.parse_known_args()
    args.conv_chr = args.conv_chr if args.conv_chr else args.annot_chr
    pretty.print_namespace(args, name_width=20); print()
    args._func(args)
