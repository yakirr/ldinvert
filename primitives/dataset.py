from __future__ import print_function, division
import json
import math
from pybedtools import BedTool
from pysnptools.snpreader import Bed
from pysnptools.util import IntRangeSet
import paths


class Dataset(object):
    def __init__(self, name, slice_size=10000):
        self.name = name
        self.__dict__.update(
                json.load(open(paths.metadata + 'datasets.json'))[name])
        bedfile = self.genotypes_bedfile()
        self.N, self.M = bedfile.iid_count, bedfile.sid_count
        self.slice_size = slice_size

    def genotypes_bedfile(self):
        return Bed(self.path + self.bfile)

    def snp_coords(self):
        return BedTool(self.path + self.bfile + '.ucscbed')

    def all_snps(self):
        return IntRangeSet((0, self.M))

    # buffer_size denotes the desired buffer to add to either side of the disjoint slices
    def slices(self, start=0, end=None, buffer_size=0):
        if not end: end = self.M
        num_slices = int(math.ceil((end - start) / self.slice_size))
        for i in xrange(num_slices):
            yield self.get_slice(i, start=start, end=end, buffer_size=buffer_size)

    def get_slice(self, slice_index, start=0, end=None, buffer_size=0):
        if not end: end = self.M
        bufferless_start = start + slice_index * self.slice_size
        bufferless_end = min(start + (slice_index+1)*self.slice_size, end)
        return (max(bufferless_start - buffer_size, start),
                min(bufferless_end + buffer_size, end))

    def get_standardized_genotypes(self, r):
        genotypes = self.genotypes_bedfile()[:,r[0]:r[1]].read()
        genotypes.standardize(); genotypes.standardize()
        return genotypes.val

if __name__ == '__main__':
    d = Dataset('tinyGERA')
    print(d.name)
    print(d.M)
    print(d.N)
    print(len(d.snp_coords()))

    d.slice_size = 300
    for s in d.slices(buffer_size=10):
        print(s)
