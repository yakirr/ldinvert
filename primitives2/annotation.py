from __future__ import print_function, division
import numpy as np
import pandas as pd
from pyutils import memo
import paths

class Annotation(object):
    def __init__(self, name, path=paths.annotations, signed=True):
        self.stem = path+name
        self.signed = signed

    @property
    def start_colindex(self):
        return 6 if self.signed else 4

    def filestem(self, chrnum=''):
        return '{}{}'.format(self.stem, chrnum)
    def filename(self, chrnum):
        if self.signed:
            return self.filestem(chrnum) + '.sannot.gz'
        else:
            return self.filestem(chrnum) + '.annot.gz'
    def sqnorm_filename(self, chrnum):
        return self.filestem(chrnum) + '.sqnorm'
    def size_filename(self, chrnum):
        return self.filestem(chrnum) + '.M'

    @memo.memoized
    def df(self, chrnum):
        return pd.read_csv(self.filename(chrnum),
                compression='gzip', header=0, sep='\t')

    def names(self, chrnum):
        return self.df(chrnum).columns.values[self.start_colindex :]

    def num_snps(self, chrnum):
        return len(self.df(chrnum))

    @memo.memoized
    def sqnorms(self, chrnum):
        return pd.read_csv(self.sqnorm_filename(chrnum), names=self.names(chrnum), sep='\t')
    @memo.memoized
    def sizes(self, chrnum):
        return pd.read_csv(self.size_filename(chrnum), names=self.names(chrnum), sep='\t')
    def total_sqnorms(self, chromosomes):
        return sum([self.sqnorms(c) for c in chromosomes])
    def total_sizes(self, chromosomes):
        return sum([self.sizes(c) for c in chromosomes])


if __name__ == '__main__':
    a = Annotation('1000G3.wim5u/mock/')
    print(a.names(22))
    print(a.df(22))
    print(a.sqnorms(22))
    print(a.sizes(22))
    print(a.total_sqnorms([1,22]))
    print(a.total_sizes([1,22]))
