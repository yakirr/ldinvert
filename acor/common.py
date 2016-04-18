from __future__ import print_function, division
import os
import pandas as pd
import numpy as np
import gzip
from pybedtools import BedTool
import primitives.dataset as prd
import primitives.genome as prg


def get_ldscores(ldscoresfile):
    print('reading ldscores')
    return pd.read_csv(ldscoresfile, header=0, sep='\t', compression='gzip')

def get_ldscores_allchr(ldscoresfile_chr, chromosomes):
    print('reading ldscores')
    ldscores = pd.DataFrame()
    for chrnum in chromosomes:
        print('\t', chrnum)
        myldscores = get_ldscores(ldscoresfile_chr + str(chrnum) + '.l2.ldscore.gz')
        ldscores = ldscores.append(myldscores)
    return ldscores

def get_sumstats(sumstatsfile):
    print('reading sumstats and filtering by sample size')
    sumstats = pd.read_csv(sumstatsfile, header=0, sep='\t',
            compression='gzip')
    N90 = np.percentile(sumstats['N'], 90)
    Nthresh = 0.9 * N90
    print('\tthreshold sample size:', Nthresh)
    print('\toriginally at', len(sumstats), 'SNPs')
    sumstats = sumstats.loc[
            sumstats['N'] >= Nthresh]
    print('\tafter filtering by N, now at', len(sumstats), 'SNPs')
    return sumstats

def get_annot(annotfilename):
    annot = pd.read_csv(annotfilename, header=0, sep='\t', compression='gzip')
    return annot, annot.columns[6:].values.tolist()

def get_annots(annotfilenames):
    print('reading annot files', annotfilenames)
    annot, annot_names = get_annot(annotfilenames[0])
    for annotfilename in annotfilenames[1:]:
        newannot, newannot_names = get_annot(annotfilename)
        toflip = annot['A1'] == newannot['A2']
        if np.sum(toflip) > 0:
            raise Exception('all annotations must have the same allele coding')
        annot[newannot_names] = newannot[newannot_names]
        annot_names += newannot_names

    print('\tannotation names:', annot_names)
    print('\tannotation contains', len(annot), 'SNPs')
    print('\tannotation supported on',
            np.sum(annot[annot_names].values != 0, axis=0), 'SNPs')
    print('\tsquared norm of annotation is',
            np.linalg.norm(annot[annot_names].values, axis=0)**2)
    return annot, annot_names

def mult_by_R_ldblocks(V, (refpanel, chrnum), ld_breakpoints, mhcpath):
    print('\tloading ld breakpoints and MHC')
    breakpoints = BedTool(ld_breakpoints)
    mhc = BedTool(mhcpath)
    print('\tconstructing SNP partition')
    blocks = prg.SnpPartition(refpanel.ucscbed(chrnum), breakpoints, mhc)

    print('\tdoing multiplication')
    result = np.zeros(V.shape)
    for r in blocks.ranges():
        print('\tXTXV', r[0], r[1], 'of', refpanel.M(chrnum))
        X = refpanel.stdX(chrnum, r)
        result[r[0]:r[1],:] = X.T.dot(X.dot(V[r[0]:r[1],:]))
    return result / refpanel.N()

def convolve(df, cols_to_convolve, (refpanel, chrnum), ld_breakpoints, mhcpath):
    print('\trefpanel contains', refpanel.M(chrnum), 'SNPs')
    print('\tmerging df and refpanel')
    if len(df) != len(refpanel.bim_df(chrnum)):
        print('df SNPs and refpanel SNPs must be the same')
        raise Exception()
    refwithdf = refpanel.bim_df(chrnum).merge(df, how='left', on=['SNP'])

    print('\tconvolving')
    RV = mult_by_R_ldblocks(refwithdf[cols_to_convolve].values, (refpanel, chrnum),
            ld_breakpoints, mhcpath)

    newnames = [n+'.conv1' for n in cols_to_convolve]
    for i, n1 in enumerate(newnames):
        df[n1] = RV[:,i]

    return df, newnames


def get_conv(convfilename, annotfilename, (refpanel, chrnum), ld_breakpoints, mhcpath):
    if not os.path.exists(convfilename):
        print('\tconv file', convfilename, 'not found. creating...')

        print('\t\treading annot', annotfilename)
        annot, annot_names = get_annot(annotfilename)

        convolved, newnames = convolve(annot, annot_names, (refpanel, chrnum),
                ld_breakpoints, mhcpath)

        def save(convolved, newnames, out):
            print('\t\tsaving')
            if out[-3:] != '.gz':
                print('\t\toutfile must end in .gz. file not saved.')
                return
            with gzip.open(out, 'w') as f:
                convolved[['CHR','BP','SNP','CM','A1','A2'] + newnames].to_csv(
                        f, index=False, sep='\t')
        save(convolved, newnames, convfilename)

    conv = pd.read_csv(convfilename, header=0, sep='\t', compression='gzip')
    return conv, conv.columns.values[6:]

def get_convs(convfilenames, annotfilenames, (refpanel, chrnum), ld_breakpoints,
        mhcpath):
    print('loading conv files', convfilenames)
    if len(convfilenames) != len(annotfilenames):
        raise Exception('\tERROR: the list of annot files and conv files must match')

    conv, conv_names = get_conv(convfilenames[0], annotfilenames[0],
            (refpanel, chrnum), ld_breakpoints, mhcpath)
    for convfilename, annotfilename in zip(convfilenames[1:], annotfilenames[1:]):
        newconv, newconv_names = get_conv(convfilename, annotfilename,
                (refpanel, chrnum), ld_breakpoints, mhcpath)
        toflip = conv['A1'] == newconv['A2']
        if np.sum(toflip) > 0:
            raise Exception('all conv files must have the same allele coding')
        conv[newconv_names] = newconv[newconv_names]
        conv_names += newconv_names

    return conv, conv_names
