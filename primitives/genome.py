from __future__ import print_function, division
import numpy as np
from glob import glob
from pybedtools import BedTool
from pysnptools.util import IntRangeSet
import paths


class GenomicSubset(object):
    def __init__(self, name, path=paths.genome_subsets):
        self.bedtool = BedTool(path + name + '.bed').sort()
        self.name = name

    @classmethod
    def whole_genome(cls, assembly='hg18'):
        return cls(assembly, path=paths.reference)

class SnpSubset(object):
    def __init__(self, genomic_subset, snp_coords):
        # use bedtools to create an indicator vector for the snps membership in the subset
        indicator = snp_coords.intersect(genomic_subset.bedtool, c=True)
        self.irs = IntRangeSet(np.flatnonzero(
            np.array([int(snp.name) for snp in indicator])))

    def num_snps(self):
        return len(self.irs)

if __name__ == '__main__':
    from dataset import Dataset
    d = Dataset('tinyGERA')
    gs = GenomicSubset('tiny')
    ss = SnpSubset(gs, d.snp_coords())
    print(ss.num_snps())
    print(GenomicSubset.whole_genome().bedtool)
