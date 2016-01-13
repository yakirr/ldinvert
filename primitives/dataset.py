from __future__ import print_function, division
import json
import numpy as np
import math
from pybedtools import BedTool
from pysnptools.snpreader import Bed
from pysnptools.util import IntRangeSet
import paths


class Dataset(object):
    def __init__(self, name, forced_M=None):
        self.name = name
        self.__dict__.update(
                json.load(open(paths.metadata + 'datasets.json'))[name])
        if forced_M is None:
            self.__bedfile = Bed(self.path + self.bfile)
        else:
            self.__bedfile = Bed(self.path + self.bfile)[:,:forced_M]
        self.N = self.genotypes_bedfile().iid_count
        self.M = self.genotypes_bedfile().sid_count

    def genotypes_bedfile(self):
        return self.__bedfile

    def snp_coords(self):
        return BedTool(self.path + self.bfile + '.ucscbed')

    def all_snps(self):
        return IntRangeSet((0, self.M))

    # buffer_size denotes the desired buffer to add to either side of the disjoint slices
    def slices(self, start=0, end=None, slice_size=10000, buffer_size=0):
        if not end: end = self.M
        num_slices = int(math.ceil((end - start) / slice_size))
        for i in xrange(num_slices):
            yield self.get_slice(
                    i, start=start, end=end, slice_size=slice_size, buffer_size=buffer_size)

    def get_slice(self, slice_index, start=0, end=None, slice_size=10000, buffer_size=0):
        if not end: end = self.M
        bufferless_start = start + slice_index * slice_size
        bufferless_end = min(start + (slice_index+1)*slice_size, end)
        return (max(bufferless_start - buffer_size, start),
                min(bufferless_end + buffer_size, end))

    def get_standardized_genotypes(self, s, indivs=None):
        if indivs is None:
            genotypes = self.genotypes_bedfile()[:,s[0]:s[1]].read()
        else:
            genotypes = self.genotypes_bedfile()[indivs, s[0]:s[1]].read()
        genotypes.standardize(); genotypes.standardize()
        return genotypes.val

    def random_indivs(self, Nref, replace=False):
        return np.random.choice(self.N, size=Nref, replace=replace)

if __name__ == '__main__':
    d = Dataset('tinyGERA')
    print(d.name)
    print(d.M)
    print(d.N)
    print(len(d.snp_coords()))

    indivs = d.random_indivs(5000)
    print(len(indivs), 'individuals subsampled')

    for s in d.slices(slice_size=300, buffer_size=10):
        print(d.get_standardized_genotypes(s, indivs=indivs).shape)
