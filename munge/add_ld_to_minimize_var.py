from __future__ import print_function, division
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--subset', type=str, required=True,
        help='the name of the genomic subset to modify')
parser.add_argument('--R2_threshold', type=float, required=True,
        help='the R2 threshold to use')
parser.add_argument('--path_to_R', type=str, required=False, default=None,
        help='path to the relevant ld matrix')
parser.add_argument('--dataset', type=str, required=True,
        help='name of the dataset corresponding to the ld matrix')
args = parser.parse_args()

import numpy as np
import pickle
from pybedtools import BedTool, Interval
from pysnptools.util import IntRangeSet
import paths
from primitives import *

def var(R, Rinv, A):
    RA = R.copy()
    A_ = np.array(list(set(range(len(R))) - set(A)))
    RA[A_] = 0
    RA.T[A_] = 0
    RinvRA = np.dot(Rinv, RA)
    # return np.trace(RinvRA.dot(RinvRA))
    return np.einsum('ij,ji', RinvRA, RinvRA)

def get_high_ld_snps(subset, R, Rinv):
    to_beat = var(R, Rinv, subset)
    add = np.argmin([var(R, Rinv, subset + [i]) for i in range(len(R))])
    if var(R, Rinv, subset + [add]) < to_beat - 10:
        print('\t\tvariance is', var(R, Rinv, subset + [add]), 'from', to_beat)
        return [add]
    else:
        print('\t\tvariance remains at', var(R, Rinv, subset))
        return []

def interval_from_range(r):
    return Interval(
            'chr'+str(int(d.genotypes_bedfile.pos[r[0]][0])),
            int(d.genotypes_bedfile.pos[r[0]][2]),
            int(d.genotypes_bedfile.pos[r[1]][2])-1)

d = Dataset(args.dataset)
A = SnpSubset(d, GenomicSubset(args.subset).bedtool)
if args.path_to_R is not None:
    R = pickle.load(open(args.path_to_R))
else:
    R = None

newA = IntRangeSet()
for r in A.expanded_by(0.003).irs.ranges():
    S = IntRangeSet([a-r[0] for a in A.irs & IntRangeSet(r)])
    print(r, 'analyzing', len(S), 'snps')
    if R is None:
        X = d.get_standardized_genotypes(r)
        cov = X.T.dot(X) / d.N
    else:
        cov = R.ranges_to_arrays[r]

    while True:
        new = get_high_ld_snps(S, cov, np.linalg.inv(cov))
        if len(new) == 0:
            break
        else:
            print('\tadding', len(new), 'snps')
            print('\t\tbefore', S)
            S += new
            print('\t\tafter', S)
    newA += IntRangeSet([s+r[0] for s in S])

b = BedTool([interval_from_range(r) for r in newA.ranges()])
print(b)
b.saveas(paths.genome_subsets + args.subset + '.R2ge{:0.2}.bed'.format(args.R2_threshold))
