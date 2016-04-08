from __future__ import print_function, division
import os
import pandas as pd
import numpy as np
import gzip
from pybedtools import BedTool
from primitives import Dataset, SnpPartition


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

def mult_by_R_ldblocks(V, refpanel, ld_breakpoints):
    print('\tloading ld breakpoints')
    breakpoints = BedTool(ld_breakpoints)
    print('\tconstructing SNP partition')
    blocks = SnpPartition(refpanel, breakpoints, remove_mhc=True)

    print('\tdoing multiplication')
    result = np.zeros(V.shape)
    for r in blocks.ranges():
        print('\tXTXV', r[0], r[1], 'of', refpanel.M)
        X = refpanel.get_standardized_genotypes(r)
        result[r[0]:r[1],:] = X.T.dot(X.dot(V[r[0]:r[1],:]))
    return result / refpanel.N

def convolve(df, cols_to_convolve, refpanel_bfile, ld_breakpoints):
    print('\treading refpanel', refpanel_bfile)
    refpanel = Dataset('ref', bfile=refpanel_bfile, reference_genome='hg19')
    print('\trefpanel contains', refpanel.M, 'SNPs')
    print('\tmerging df and refpanel')
    if len(df) != len(refpanel.bim):
        print('df SNPs and refpanel SNPs must be the same')
        raise Exception()
    refwithdf = refpanel.bim.merge(df, how='left', on=['SNP'])

    print('\tconvolving')
    RV = mult_by_R_ldblocks(refwithdf[cols_to_convolve].values, refpanel, ld_breakpoints)

    newnames = [n+'.conv1' for n in cols_to_convolve]
    for i, n1 in enumerate(newnames):
        df[n1] = RV[:,i]

    return df, newnames


def get_conv(convfilename, annotfilename, refpanel_bfile, ld_breakpoints):
    if not os.path.exists(convfilename):
        print('\tconv file', convfilename, 'not found. creating...')

        print('\t\treading annot', annotfilename)
        annot, annot_names = get_annot(annotfilename)

        convolved, newnames = convolve(annot, annot_names, refpanel_bfile, ld_breakpoints)
        import pdb; pdb.set_trace()

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

def get_convs(convfilenames, annotfilenames, refpanel_bfile, ld_breakpoints):
    print('loading conv files', convfilenames)
    if len(convfilenames) != len(annotfilenames):
        raise Exception('\tERROR: the list of annot files and conv files must match')

    conv, conv_names = get_conv(convfilenames[0], annotfilenames[0],
            refpanel_bfile, ld_breakpoints)
    for convfilename, annotfilename in zip(convfilenames[1:], annotfilenames[1:]):
        newconv, newconv_names = get_conv(convfilename, annotfilename,
                refpanel_bfile, ld_breakpoints)
        toflip = conv['A1'] == newconv['A2']
        if np.sum(toflip) > 0:
            raise Exception('all conv files must have the same allele coding')
        conv[newconv_names] = newconv[newconv_names]
        conv_names += newconv_names

    return conv, conv_names









def results_filename(annot_stems, sumstats_stem, chrnum=None):
    annot_name = '-'.join([annot_stem.split('/')[-1] for annot_stem in annot_stems])
    pheno_name = sumstats_stem.split('/')[-1]
    stem = '/groups/price/yakir/results/2016.03.29_signed/' + \
            pheno_name + '.' + annot_name
    if chrnum is None:
        return stem
    else:
        return stem + '.' + str(chrnum)

# def get_annot(args):
#     print('reading annot files')
#     annot = None
#     for annot_stem in args.annot_stems:
#         print('\t', annot_stem)
#         if annot is None:
#             annot = pd.read_csv('{}.{}.annot.gz'.format(annot_stem, args.chrnum),
#                     header=0, sep='\t', compression='gzip')
#         else:
#             newannot = pd.read_csv('{}.{}.annot.gz'.format(annot_stem, args.chrnum),
#                     header=0, sep='\t', compression='gzip')
#             newannot_names = newannot.columns[4:].values.tolist()
#             newannot = newannot[['SNP'] + newannot_names]
#             annot = annot.merge(newannot, how='left', on='SNP')
#     annot_names = annot.columns[4:]
#     print('\tannotation names:', annot_names.values)
#     print('\tannotations supported on',
#             np.sum(annot[annot_names].values != 0, axis=0), 'SNPs')
#     return annot, annot_names
#
# def get_refpanel(args):
#     print('reading refpanel and refpanel bim')
#     refpanel = Dataset(args.refpanel + '.' + str(args.chrnum))
#     refpanel_bim = pd.read_csv(refpanel.genotypes_bedfile.filename + '.bim',
#             sep='\t',
#             names=['CHR', 'SNP', 'cM', 'BP', 'A1', 'A2'])
#     refpanel_bim['refpanelINDEX'] = np.arange(len(refpanel_bim))
#     print('\trefpanel contains', len(refpanel_bim), 'SNPs')
#     return refpanel, refpanel_bim
#
# def merge_files(refpanel_bim, sumstats, annot, annot_names, ldscores=None,
#         get_full_annot=False):
#     print('merging files')
#     print('\tmerging refpanel and sumstats')
#     data = refpanel_bim.merge(sumstats, how='left', on=['SNP'])
#     print('\tmerging annot')
#     data = data.merge(annot, how='left', on=['SNP'])
#     if get_full_annot:
#         V = data[annot_names].values
#
#     data = data.loc[
#             data['N'].notnull()]
#     print('\tnumber of SNPs in common between refpanel and sumstats is', len(data))
#     print('\tannotations now supported on',
#             np.sum(data[annot_names].values != 0, axis=0), 'SNPs')
#
#     if ldscores is not None:
#         print('\tmerging ldscores')
#         data = data.merge(ldscores, how='left', on=['SNP'])
#
#     print('\tflipping sumstats to line up with refpanel')
#     flip = data['A1_x'] == data['A2_y']
#     data.ix[flip, 'Z'] *= -1
#
#     print('\tcreating numpy arrays')
#     N = np.min(data['N'].values)
#     alphahat = data['Z'].values / np.sqrt(data['N'].values)
#     l2_ss = None if ldscores is None else data['L2'].values
#     if not get_full_annot:
#         V = data[annot_names].values
#     refpanel_indices = data['refpanelINDEX'].values
#     print('\tN=', N, 'alphahat.shape=', alphahat.shape, 'V.shape=', V.shape,
#             'refpanel_indices.shape=', refpanel_indices.shape)
#
#     return N, alphahat, V, refpanel_indices, l2_ss
#
# def mult_by_R_ldblocks(V, refpanel, refpanel_indices=None):
#     breakpoints = BedTool(paths.reference + 'pickrell_breakpoints.hg19.eur.bed')
#     blocks = SnpPartition(refpanel, breakpoints, remove_mhc=True)
#
#     if refpanel_indices is not None:
#         Vwithzeros = np.zeros((refpanel.M, V.shape[1]))
#         Vwithzeros[refpanel_indices] = V
#     else:
#         Vwithzeros = V
#     result = np.zeros(Vwithzeros.shape)
#     for r in blocks.ranges():
#         print('\tXTXV', r[0], r[1], 'of', refpanel.M)
#         X = refpanel.get_standardized_genotypes(r)
#         result[r[0]:r[1],:] = X.T.dot(X.dot(Vwithzeros[r[0]:r[1],:]))
#     if refpanel_indices is not None:
#         return result[refpanel_indices] / refpanel.N
#     else:
#         return result / refpanel.N
