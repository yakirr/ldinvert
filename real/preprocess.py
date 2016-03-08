from __future__ import print_function, division
import argparse
import numpy as np
from methods import ldinvert2
from primitives import Dataset
import paths
from pyutils import bsub

def main(args):
    est = ldinvert2.Ldinvert(refpanel='UK10Khg19.'+str(args.chrom),
            breakpointsfile='pickrell_breakpoints.hg19.eur.bed',
            region='50',
            Lambda=0.03,
            pop_size=10**15,
            num_chunks=5,
            chunk=args.chunk)

    est.preprocess()

def submit(args):
    my_args = ['main',
            '--chrom', '$LSB_JOBINDEX']
    d = Dataset('UK10Khg19.22')
    outfilepath = d.auxfiles_path + '../' + \
            '%I/.preprocess.out'
    bsub.submit(
            ['python', '-u', paths.code + 'real/preprocess.py'] + my_args,
            outfilepath,
            jobname='preprocess[1-22]',
            memory_GB=16)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparser_main, subparser_submit = bsub.add_main_and_submit(parser,
            main,
            submit)
    subparser_main.add_argument('--chrom', type=int, required=True,
            help='the chromosome to preprocess')
    subparser_main.add_argument('--chunk', type=int, required=True,
            help='the chunk of the chromosome to process (1-indexed)')

    bsub.choose_parser_and_run(parser)
