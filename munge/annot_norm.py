from __future__ import print_function, division
import numpy as np
import pandas as pd
import paths
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--annot_name', required=True, type=str,
        help='without the .annot.gz suffix')

args = parser.parse_args()

a = pd.read_csv(args.annot_name+'.annot.gz', header=0, sep='\t', compression='gzip')
names = a.columns[4:].values
result = np.linalg.norm(a[names], axis=0)**2
print(names)
print(result)

with open(args.annot_name + '.norm', 'w') as f:
    print('\t'.join(names), file=f)
    print('\t'.join([str(r) for r in result]), file=f)
