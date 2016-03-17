from __future__ import print_function, division
import argparse
import pandas as pd
import numpy as np
import gzip
from primitives import Dataset, SnpSubset, GenomicSubset
from pyutils import bsub
import paths


def main(args):
    refpanel = Dataset(args.refpanel + '.' + str(args.chrnum))
    annot_filename = '{}.{}.annot.gz'.format(args.annot_stem, args.chrnum)
    result_filename = '{}.maf5p.{}.annot.gz'.format(args.annot_stem, args.chrnum)

    print('reading annot')
    annot = pd.read_csv(annot_filename, compression='gzip', sep='\t', header=0)
    name = annot.columns[-1]
    names = annot.columns.values

    print('reading frq')
    refpanel_frq = pd.read_csv(refpanel.genotypes_bedfile.filename + '.frq',
            delim_whitespace=True,
            header=0)
    refpanel_frq = refpanel_frq[['SNP', 'MAF']]

    print('merging')
    annot = annot.merge(refpanel_frq, how='left', on=['SNP'])

    print('before filtering, |A| =', np.sum(annot[name]))
    annot.ix[annot['MAF'] > 0.05, name] = 0
    print('after filtering, |A| =', np.sum(annot[name]))

    # write output
    print('writing output')
    with gzip.open(result_filename, 'wt') as f:
        annot[names].to_csv(f, index=False, sep='\t')


def submit(args):
    my_args = ['--annot_stem', args.annot_stem,
            '--refpanel', args.refpanel,
            'main',
            '--chrnum', '$LSB_JOBINDEX']
    outfilepath = \
            args.annot_stem + '.maf5p.%I.creation.log'

    bsub.submit(
            ['python', '-u', paths.code + 'real/select_rare_only.py'] + my_args,
            outfilepath,
            jobname='selectrare[1-22]')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--annot_stem', type=str, required=True,
            help='path to annot file, not including chromosome number and extension.  ' + \
                    'For example: path/to/annot')
    parser.add_argument('--refpanel', type=str, required=True, help='reference panel to use')

    main_parser, submit_parser = bsub.add_main_and_submit(parser, main, submit)
    main_parser.add_argument('--chrnum', type=int, required=True, help='chromosome to analyze')

    bsub.choose_parser_and_run(parser)
