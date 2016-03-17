from __future__ import print_function, division
import argparse
import pandas as pd
import numpy as np
import gzip
from primitives import Dataset, SnpSubset, GenomicSubset
from pyutils import bsub
import paths


def main(args):
    annot_filename = '{}.{}.annot.gz'.format(args.annot_stem, args.chrnum)
    result_filename = '{}.random.{}.annot.gz'.format(args.annot_stem, args.chrnum)

    print('reading annot')
    annot = pd.read_csv(annot_filename, compression='gzip', sep='\t', header=0)
    name = annot.columns[-1]

    random_signs = (-1)**np.random.randint(2, size=len(annot))
    annot[name] *= random_signs

    # write output
    print('writing output')
    with gzip.open(result_filename, 'wt') as f:
        annot.to_csv(f, index=False, sep='\t')


def submit(args):
    my_args = ['--annot_stem', args.annot_stem,
            'main',
            '--chrnum', '$LSB_JOBINDEX']
    outfilepath = \
            args.annot_stem + '.random.%I.creation.log'

    bsub.submit(
            ['python', '-u', paths.code + 'munge/randomize_annot_signs.py'] + my_args,
            outfilepath,
            jobname='randomize[1-22]', debug=args.debug)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--annot_stem', type=str, required=True,
            help='path to annot file, not including chromosome number and extension.  ' + \
                    'For example: path/to/annot')

    main_parser, submit_parser = bsub.add_main_and_submit(parser, main, submit)
    main_parser.add_argument('--chrnum', type=int, required=True, help='chromosome to analyze')
    submit_parser.add_argument('-debug', action='store_true', default=False,
            help='do not actually submit the jobs')

    bsub.choose_parser_and_run(parser)
