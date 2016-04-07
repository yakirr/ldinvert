from __future__ import print_function, division
import gzip
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
from pybedtools import BedTool
from primitives import Dataset, SnpPartition
import common


def get_annot(annotfilename):
    print('reading annot file', annotfilename)
    annot = pd.read_csv(annotfilename, header=0, sep='\t', compression='gzip')
    annot_names = annot.columns[4:].values.tolist()
    print('\tannotation names:', annot_names)
    print('\tannotation contains', len(annot), 'SNPs')
    print('\tannotation supported on',
            np.sum(annot[annot_names].values != 0, axis=0), 'SNPs')
    print('\tsquared norm of annotation is',
            np.linalg.norm(annot[annot_names].values, axis=0)**2)
    return annot, annot_names

def get_refpanel(bfile):
    print('reading refpanel', bfile)
    refpanel = Dataset('ref', bfile=bfile, reference_genome='hg19')
    # refpanel_bim['refpanelINDEX'] = np.arange(len(refpanel_bim))
    print('\trefpanel contains', refpanel.M, 'SNPs')
    return refpanel

def merge_files(annot, refpanel_bim, ldscores):
    print('merging files')
    print('\tmerging annot')
    if len(annot) != len(refpanel_bim):
        print('annotation SNPs and refpanel SNPs must be the same')
        raise Exception()
    data = refpanel_bim.merge(annot, how='left', on=['SNP'])

    print('\tmerging in ld scores')
    if len(ldscores) != len(refpanel_bim):
        print('ldscores SNPs and refpanel SNPs must be the same')
        raise Exception()
    data = data.merge(ldscores, how='left', on=['SNP'])

    return data

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

def save(RV, R3V, annot, annot_names, refpanel, out):
    print('saving')
    if out[-3:] != '.gz':
        print('\toutfile must end in .gz. file not saved.')
        return
    newnames = [[n+'.conv',n+'.conv3'] for n in annot_names]
    newnames_flat = sum(newnames, [])

    for i, [n1, n3] in enumerate(newnames):
        refpanel.bim[n1] = RV[:,i]
        refpanel.bim[n3] = R3V[:,i]
    with gzip.open(out, 'w') as f:
        refpanel.bim[['CHR','BP','SNP'] + newnames_flat].to_csv(f, index=False, sep='\t')

def conv(annot=None, bfile=None, ld_breakpoints=None, out=None, ldscores=None):
    annot, annot_names = get_annot(annot)
    refpanel = get_refpanel(bfile)
    data = merge_files(annot, refpanel.bim, common.get_ldscores(ldscores))

    print('computing RV')
    RV = mult_by_R_ldblocks(data[annot_names].values, refpanel, ld_breakpoints)

    print('computing R2V')
    R2V = mult_by_R_ldblocks(RV, refpanel, ld_breakpoints)

    print('computing R3V')
    R3V = mult_by_R_ldblocks(R2V, refpanel, ld_breakpoints)

    save(RV, R3V, annot, annot_names, refpanel, out)

    return RV, R3V
