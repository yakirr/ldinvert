from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('--bfile-chr')
parser.add_argument('--chrom-list', nargs='+',
        default=map(str,range(1,23)))
parser.add_argument('--p1', type = float)
parser.add_argument('--v', type = float)
parser.add_argument('--out-prefix', default='')
args = parser.parse_args()

np.random.seed(0)

for chrom in args.chrom_list:
    print(chrom)
    bfile = args.bfile_chr + chrom
    df = pd.read_csv(bfile + '.bim', delim_whitespace = True, header = None, usecols = [0,1,2,3,4,5], names = ['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
    a = np.zeros(len(df))
    mask = np.random.binomial(1, args.p1, len(df)).astype(bool)
    a[mask] = np.random.normal(0, args.v**0.5, sum(mask))
    df['ANNOT'] = a
    df = df[['CHR', 'BP', 'SNP', 'CM', 'A1', 'A2', 'ANNOT']]

    with gzip.open(args.out_prefix + chrom + '.sannot.gz', 'w') as f:
        df.to_csv(f, sep = "\t", index = False)
