from __future__ import print_function, division
import copy
import numpy as np
from glob import glob
from pybedtools import BedTool
from pysnptools.util import IntRangeSet
import paths


class GenomicSubset(object):
    def __init__(self, name, path=paths.genome_subsets, assembly='hg18'):
        self.assembly = assembly
        self.name = name
        self.bedtool = BedTool(path + name + '.bed').sort()

        # Intersect the pathway with the appropriate genome build
        # TODO: this step should be unnecessary if the pathways are correct
        if name != self.assembly:
            self.bedtool = GenomicSubset.reference_genome(
                    self.assembly).bedtool.intersect(self.bedtool).sort().saveas()

    def expand_by(self, expansion_in_each_direction_Mb):
        window_size_str = str(expansion_in_each_direction_Mb) + 'Mb'
        print('total size before window addition:', self.bedtool.total_coverage(), 'bp')

        # compute the flanks
        # TODO: use 1cM instead of 1Mb
        print('computing flanks')
        flanks = self.bedtool.flank(
            genome=self.assembly,
            b=expansion_in_each_direction_Mb*1000000).sort().merge().saveas()

        # compute the union of the flanks and the pathway
        print('computing union')
        union = self.bedtool.cat(flanks, postmerge=False).sort()
        merged = union.merge().saveas()
        print('total size after window addition:', merged.total_coverage(), 'bp')
        self.bedtool = merged

    @classmethod
    def reference_genome(cls, assembly='hg18'):
        return GenomicSubset(assembly, path=paths.reference, assembly=assembly)

    @classmethod
    def whole_genome(cls, assembly='hg18'):
        return cls(assembly, path=paths.reference)

class SnpSubset(object):
    def __init__(self, genomic_subset, dataset):
        # use bedtools to create an indicator vector for the snps membership in the subset
        self.dataset = dataset
        indicator = dataset.snp_coords().intersect(genomic_subset.bedtool, c=True)
        self.irs = IntRangeSet(np.flatnonzero(
            np.array([int(snp.name) for snp in indicator])))

    def num_snps(self):
        return len(self.irs)

    def expand_by(self, expansion_in_each_direction, units='Morgans'):
        result = IntRangeSet()
        for r in self.irs.ranges():
            if units == 'Mb':
                result += (max(0, r[0] - int(10**6 * expansion_in_each_direction)),
                    r[1] + int(10**6 * expansion_in_each_direction))
            elif units == 'Morgans':
                result += self.dataset.buffer_around_slice(r, expansion_in_each_direction)
        self.irs = result

    def expanded_by(self, expansion_in_each_direction, units='Morgans'):
        result = copy.copy(self)
        result.expand_by(expansion_in_each_direction, units=units)
        return result

if __name__ == '__main__':
    from dataset import Dataset
    d = Dataset('tinyGERA')
    gs = GenomicSubset('tiny')
    ss = SnpSubset(gs, d)
    print(ss.num_snps())
    print(GenomicSubset.whole_genome().bedtool)
