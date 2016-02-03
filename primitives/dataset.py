from __future__ import print_function, division
import json
import os
import numpy as np
import pandas as pd
import bisect
import math
from pybedtools import BedTool
from pysnptools.snpreader import Bed
from pysnptools.util import IntRangeSet
import paths


class Dataset(object):
    def __init__(self, name, forced_M=None, chrnum=None):
        self.name = name
        self.__dict__.update(
                json.load(open(paths.metadata + 'datasets.json'))[name])
        if chrnum:
            self.name += '|chr' + str(int(chrnum))
            self.bfile += '.' + str(int(chrnum))

        if forced_M is None:
            self.__bedfile = Bed(self.path + self.bfile)
        else:
            self.__bedfile = Bed(self.path + self.bfile)[:,:forced_M]
        self.N = self.genotypes_bedfile().iid_count
        self.M = self.genotypes_bedfile().sid_count
        # chrom 1 starts at chrom_boundaries[0] and ends before chrom_boundaries[1])
        self.chrom_boundaries = np.searchsorted(self.genotypes_bedfile().pos[:,0],
                np.arange(23) +1 )

    def genotypes_bedfile(self):
        return self.__bedfile

    def snp_coords(self):
        ucscbedfilename = self.path + self.bfile + '.ucscbed'
        if os.path.exists(ucscbedfilename):
            return BedTool(self.path + self.bfile + '.ucscbed')
        else:
            print('warning: ucscbedfile not found:', ucscbedfilename)
            return None

    def all_snps(self):
        return IntRangeSet((0, self.M))

    def chromosomes(self):
        return np.unique(self.genotypes_bedfile().pos[:,0]).astype(np.int16)

    def snp_at_distance(self, snp_index, distance_in_morg, break_ties_to_right=True,
            units='Morgans'):
        if units == 'Morgans':
            i = 1
        elif units == 'bp':
            i = 2
        else:
            print('ERROR: units must be either Morgans or bp')
            return None

        snp_index = min(snp_index, self.M-1)
        coords = self.genotypes_bedfile().pos[:,i]
        chrom = self.genotypes_bedfile().pos[snp_index,0]
        chrom_coords = coords[self.chrom_boundaries[chrom-1]:self.chrom_boundaries[chrom]]

        if break_ties_to_right:
            return bisect.bisect_right(chrom_coords,
                    coords[snp_index] + distance_in_morg) - 1 + self.chrom_boundaries[chrom-1]
        else:
            return bisect.bisect_left(chrom_coords,
                    coords[snp_index] + distance_in_morg) + self.chrom_boundaries[chrom-1]

    # buffer_size denotes the desired buffer to add to either side of the disjoint slices
    def slices(self, start=0, end=None, slice_size=10000):
        if not end: end = self.M
        num_slices = int(math.ceil((end - start) / slice_size))
        for i in xrange(num_slices):
            yield self.get_slice(
                    i, start=start, end=end, slice_size=slice_size)

    def get_slice(self, slice_index, start=0, end=None, slice_size=10000):
        if not end: end = self.M
        candidate_start = start + slice_index * slice_size
        candidate_end = min(start + (slice_index+1)*slice_size, end)
        return (max(candidate_start, start),
                min(candidate_end, end))

    def buffer_around_slice(self, s, buffer_size_morg, start=0, end=None, units='Morgans'):
        if not end: end = self.M
        buffered_start = self.snp_at_distance(s[0], -buffer_size_morg,
                break_ties_to_right=False, units=units)
        buffered_end = self.snp_at_distance(s[1]-1, buffer_size_morg,
                break_ties_to_right=True, units=units) + 1
        return (max(buffered_start, start),
                min(buffered_end, end))

    def buffer_around_snp(self, snp_index, buffer_size_morg, start=0, end=None,
            units='Morgans'):
        return self.buffer_around_slice((snp_index, snp_index+1),
                buffer_size_morg,
                start=start, end=end,
                units=units)

    def get_standardized_genotypes(self, s, indivs=None):
        if indivs is None:
            genotypes = self.genotypes_bedfile()[:,s[0]:s[1]].read()
        else:
            genotypes = self.genotypes_bedfile()[indivs, s[0]:s[1]].read()
        genotypes.standardize(); genotypes.standardize()
        return genotypes.val

    def random_indivs(self, Nref, replace=False):
        return np.random.choice(self.N, size=Nref, replace=replace)

    def random_indivs_df(self, Nref, replace=False):
        indivs = self.random_indivs(Nref, replace=replace)
        all_indivs = pd.read_csv(self.genotypes_bedfile().filename + '.fam',
                delim_whitespace=True,
                header=None,
                usecols=[0,1])
        return all_indivs.ix[indivs]


if __name__ == '__main__':
    d = Dataset('tinyGERA')
    print(d.name)
    print(d.M)
    print(d.N)
    print(len(d.snp_coords()))

    indivs = d.random_indivs(5000)
    print(len(indivs), 'individuals subsampled')

    for s in d.slices(slice_size=300):
        print('slice', s)
        bs = d.buffer_around_slice(s, 10000, units='bp')
        print('with buffer:', bs)
        print('shape:', d.get_standardized_genotypes(bs, indivs=indivs).shape)
        print(d.genotypes_bedfile().pos[bs[0],2], d.genotypes_bedfile().pos[s[0],2])
