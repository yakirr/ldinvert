from __future__ import print_function, division
import argparse
import pandas as pd
import numpy as np
import gzip
from primitives import Dataset, SnpSubset, GenomicSubset
from pyutils import bsub
import paths


def main(args):
    print('reading sumstats')
    sumstats = pd.read_csv(args.sumstats_stem + '.sumstats.gz', header=0, sep='\t',
            compression='gzip')

    point_estimate = 0
    variance = 0
    for chrnum in range(1,23):
        print(chrnum, 'reading cannot file')
        cannot = pd.read_csv('{}.{}.cannot.gz'.format(args.annot_stem, chrnum),
                header=0, sep='\t', compression='gzip')
        conv_annot_name = [s for s in cannot.columns if '.CONV' in s][0]
        annot_name = conv_annot_name[:-5]

        print('\treading refpanel bim')
        refpanel = Dataset(args.refpanel + '.' + str(chrnum))
        refpanel_bim = pd.read_csv(refpanel.genotypes_bedfile.filename + '.bim',
                sep='\t',
                names=['CHR', 'SNP', 'cM', 'BP', 'A1', 'A2'])
        data = refpanel_bim.merge(sumstats, how='left', on=['SNP'])
        data = data.loc[
                data['N'].notnull()]
        data = data.merge(cannot, how='left', on=['SNP'])

        print('\tflipping sumstats to line up with refpanel')
        flip = data['A1_x'] == data['A2_y']
        data.ix[flip, 'Z'] *= -1

        N = np.mean(data.ix[:,'N'].values) #TODO: make this a threshold instead
        alphahat = data.ix[:,'Z'].values / np.sqrt(N)
        Rv = data.ix[:,conv_annot_name].values
        v = data.ix[:,annot_name].values
        with open('{}.{}.cannot.norm'.format(args.annot_stem, chrnum), 'r') as f:
            vTRv = float(f.readline())
        print(alphahat.shape)
        print(vTRv)

        my_point_estimate = alphahat.dot(v)
        my_variance = vTRv / N + point_estimate**2 / N
        print('chr', chrnum, ':', my_point_estimate, '(', np.sqrt(my_variance), ')')
        point_estimate += my_point_estimate
        variance += my_variance
    print('total :', point_estimate, '(', np.sqrt(variance), ')')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--annot_stem', type=str, required=True,
            help='path to (c)annot file, not including chromosome number and extension.  ' + \
                    'For example: path/to/annot')
    parser.add_argument('--sumstats_stem', type=str, required=True,
            help='path to sumstats file, not including .sumstats.gz extension.')
    parser.add_argument('--refpanel', type=str, required=True,
            help='the name of the reference panel, for synchronizing allele coding')

    args = parser.parse_args()
    main(args)

