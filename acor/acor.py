from __future__ import print_function, division
import argparse
import pyutils.pretty as pretty
import conv


def cor_fe_wrap(args):
    print('cor_fe')
    #TODO
    pass

def cor_re_wrap(args):
    print('cor_re')
    #TODO
    pass

def conv_wrap(args):
    print('\nCONV')
    if args.out is None:
        if '.annot' in args.annot:
            args.out = args.annot[:args.annot.rfind('.annot')] + '.conv.gz'
        else:
            args.out = args.annot + '.conv.gz'
        print('--out unspecified. Defaulting to', args.out)
    print()

    d = vars(args)
    del d['_func']
    conv.conv(**d)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    subparser_conv = subparsers.add_parser('conv')
    subparser_conv.set_defaults(_func=conv_wrap)
    subparser_conv.add_argument('--annot', required=True,
            help='path to gzipped annot file, including chromosome number if any, and suffix')
    subparser_conv.add_argument('--bfile', required=True,
            help='path to plink bfile of reference panel to use')
    subparser_conv.add_argument('--ld-breakpoints', required=True,
            help='path to UCSC bed file containing one zero-length bed interval per LD' + \
                    ' breakpoint')
    subparser_conv.add_argument('--ldscores', required=True,
            help='path to a .l2.ldscore.gz file containing a column named L2 with ld scores')
    subparser_conv.add_argument('--out', default=None,
            help='stem of output files, including chromosome number if any, but no extension')

    subparser_cor = subparsers.add_parser('cor')
    subparser_cor.add_argument('--ref-ld-chr', required=True,
            help='path to reference panel l2.ldscore.gz/conv.gz files, not including' + \
                    ' chromosome number')
    subparser_cor.add_argument('--sumstats', required=True,
            help='path to sumstats.gz file, including extension')
    subparser_cor.add_argument('--out', required=True,
            help='path to output file')
    #TODO: add ability to input separate ld scores for weights, i.e., ld-scores TO only the
    #   regression snps.
    subparser_cor_subs = subparser_cor.add_subparsers()
    subparser_cor_fe = subparser_cor_subs.add_parser('fe')
    subparser_cor_fe.set_defaults(_func=cor_fe_wrap)
    subparser_cor_re = subparser_cor_subs.add_parser('re')
    subparser_cor_re.set_defaults(_func=cor_re_wrap)

    args, _ = parser.parse_known_args()
    pretty.print_namespace(args); print()
    args._func(args)
