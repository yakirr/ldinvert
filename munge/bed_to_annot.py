from __future__ import print_function, division
import argparse
from primitives import Dataset, SnpSubset, GenomicSubset
import paths


def create_annot(args):
    path = '/'.join(args.bedfile.split('/')[:-1]) + '/'
    filename = args.bedfile.split('/')[-1]
    if filename[-4:] == '.bed':
        name = filename[:-4]
    else:
        name = filename

    gs = GenomicSubset(name, path=path)
    for chrnum in range(1,23)[::-1]:
        print('creating annot file for chr', chrnum)
        d = Dataset(args.refpanel + '.' + str(chrnum))
        sss = [SnpSubset(d, gs.restricted_to_chrom_bedtool(chrnum))]
        SnpSubset.print_subsets('{}{}.{}.annot.gz'.format(path, name, chrnum),
                sss, [name])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedfile', type=str, required=True, help='path to bed file')
    parser.add_argument('--refpanel', type=str, required=True, help='name of refpanel')

    args = parser.parse_args()
    create_annot(args)

