from __future__ import print_function, division
import numpy as np
import pandas as pd
import paths
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--annot-name', required=True, type=str,
        help='without the .sannot.gz prefix')

args = parser.parse_args()

a = pd.read_csv(args.annot_name+'.sannot.gz', header=0, sep='\t', compression='gzip')
names = a.columns[6:].values
norms = np.linalg.norm(a[names], axis=0)**2
sizes = np.sum(a[names].values !=0, axis=0)
print(names)
print(norms)
print(sizes)

with open(args.annot_name + '.sqnorm', 'w') as f:
    print('\t'.join([str(r) for r in norms]), file=f)

with open(args.annot_name + '.M', 'w') as f:
    print('\t'.join([str(r) for r in sizes]), file=f)
