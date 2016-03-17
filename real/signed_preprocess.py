from __future__ import print_function, division
import argparse
import pandas as pd
import numpy as np
import gzip
from primitives import Dataset, SnpSubset, GenomicSubset
from pyutils import bsub
import paths


def main(args):
    d = Dataset(args.refpanel + '.' + str(args.chrnum))
    annot_filename = '{}.{}.annot.gz'.format(args.annot_stem, args.chrnum)
    cannot_filename = '{}.{}.cannot.gz'.format(args.annot_stem, args.chrnum)
    cannot_norm_filename = '{}.{}.cannot.norm'.format(args.annot_stem, args.chrnum)

    annot = pd.read_csv(annot_filename, compression='gzip', sep='\t', header=0)
    name = annot.columns[-1]
    v = annot.ix[:,name].values

    #TODO: use ld blocks, possibly just those that have non-trivial intersection with the
    # nonzero entries of v
    print('computing Xv')
    Xv = np.zeros(d.N)
    for s in d.slices():
        print(s)
        X = d.get_standardized_genotypes(s)
        Xv += X.dot(v[s[0]:s[1]])

    print('computing XTXv')
    XTXv = np.zeros(d.M)
    for s in d.slices():
        print(s)
        X = d.get_standardized_genotypes(s)
        XTXv[s[0]:s[1]] = X.T.dot(Xv)

    print('computing V^TRv')
    Rv = XTXv / d.N
    vTRv = v.dot(Rv)

    # write output
    print('writing output')
    annot[name+'.CONV'] = Rv
    with gzip.open(cannot_filename, 'wt') as f:
        annot.to_csv(f, index=False, sep='\t')

    with open(cannot_norm_filename, 'w') as f:
        f.write(str(vTRv))

def submit(args):
    my_args = ['--annot_stem', args.annot_stem,
            '--refpanel', args.refpanel,
            'main',
            '--chrnum', '$LSB_JOBINDEX']
    outfilepath = \
            args.annot_stem + '.%I.cannot.log'

    bsub.submit(
            ['python', '-u', paths.code + 'real/signed_preprocess.py'] + my_args,
            outfilepath,
            jobname='preprocess[1-22]')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--annot_stem', type=str, required=True,
            help='path to annot file, not including chromosome number and extension.  ' + \
                    'For example: path/to/annot')
    parser.add_argument('--refpanel', type=str, required=True, help='reference panel to use')

    main_parser, submit_parser = bsub.add_main_and_submit(parser, main, submit)
    main_parser.add_argument('--chrnum', type=int, required=True, help='chromosome to analyze')

    bsub.choose_parser_and_run(parser)
