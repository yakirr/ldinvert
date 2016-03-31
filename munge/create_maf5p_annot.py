from __future__ import print_function, division
import argparse
import pandas as pd
import numpy as np
import gzip
from primitives import Dataset, SnpSubset, GenomicSubset
from pyutils import bsub
import paths


def main(args):
    result_filename = '{}.{}.annot.gz'.format(args.annot_stem, args.chrnum)
    name = args.annot_stem.split('/')[-1]

    print('reading refpanel bim')
    refpanel = Dataset(args.refpanel + '.' + str(args.chrnum))
    refpanel_bim = pd.read_csv(refpanel.genotypes_bedfile.filename + '.bim',
            sep='\t',
            names=['CHR', 'SNP', 'cM', 'BP', 'A1', 'A2'])
    print('\tthere are', len(refpanel_bim), 'SNPs')

    print('reading frq')
    refpanel_frq = pd.read_csv(refpanel.genotypes_bedfile.filename + '.frq',
            delim_whitespace=True,
            header=0)
    refpanel_frq = refpanel_frq[['SNP', 'MAF']]

    print('merging')
    annot = refpanel_bim.merge(refpanel_frq, how='left', on=['SNP'])
    annot[name] = 0

    print('before filtering, |A| =', np.sum(annot[name]))
    to_include = np.flatnonzero((annot['MAF'] <= 0.01).values)
    # sizes = [
    #         27285,
    #         26794,
    #         23253,
    #         18865,
    #         20449,
    #         17280,
    #         15776,
    #         13198,
    #         14079,
    #         15070,
    #         14247,
    #         12971,
    #         9039,
    #         11219,
    #         8688,
    #         11506,
    #         10187,
    #         7603,
    #         7804,
    #         6765,
    #         3747,
    #         4021
    #         ]
    # to_include = to_include[np.random.choice(len(to_include), replace=False,
        # size=sizes[args.chrnum-1])]
    annot.ix[to_include, name] = 1
    print('after filtering, |A| =', np.sum(annot[name]))
    annot.rename(columns={'cM':'CM'}, inplace = True)
    names = ['CHR', 'BP', 'SNP', 'CM', name]

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
            args.annot_stem + '%I.creation.log'

    bsub.submit(
            ['python', '-u', __file__] + my_args,
            outfilepath,
            jobname='maf1p[1-22]', debug=args.debug)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--annot_stem', type=str, required=True,
            help='path to output annot file, not including chromosome number and ' + \
                    'extension. For example: path/to/annot')
    parser.add_argument('--refpanel', type=str, required=True, help='reference panel to use')

    main_parser, submit_parser = bsub.add_main_and_submit(parser, main, submit)
    main_parser.add_argument('--chrnum', type=int, required=True, help='chromosome to create')
    submit_parser.add_argument('-debug', action='store_true', default=False)

    bsub.choose_parser_and_run(parser)
