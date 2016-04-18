from __future__ import print_function, division
import numpy as np
from pysnptools.util import IntRangeSet


class SnpSubset(object):
    def __init__(self, snps_ucscbed, ucscbed):
        # use bedtools to create an indicator vector for the snps membership in the subset
        indicator = snps_ucscbed.intersect(ucscbed, c=True)
        self.irs = IntRangeSet(np.flatnonzero(
            np.array([int(snp.name) for snp in indicator])))

    def num_snps(self):
        return len(self.irs)

class SnpPartition(object):
    def __init__(self, snps_ucscbed, breakpoints_ucscbed, mhc_ucscbed=None):
        # use bedtools to find the snp indices of the breakpoints
        self.last_snps = np.array([])
        c = snps_ucscbed.closest(
                breakpoints_ucscbed.sort(), iu=True, D='ref').saveas()
        closest_breakpoints = np.array([' '.join(i[3:-1]) for i in c] + [''])
        indices = closest_breakpoints[:-1] != closest_breakpoints[1:]
        self.last_snps = np.flatnonzero(indices)

        # store MHC info if supplied
        self.mhc_ucscbed = mhc_ucscbed
        if self.mhc_ucscbed is not None:
            self.mhc = SnpSubset(snps_ucscbed, self.mhc_ucscbed)

    def ranges(self):
        ranges = zip(np.concatenate([[0], self.last_snps]), self.last_snps+1)
        if self.mhc_ucscbed is None:
            return ranges
        else:
            return [r for r in ranges if (IntRangeSet(r) & self.mhc.irs).isempty]

    def ind_ranges_overlapping(self, irs):
        ranges = self.ranges()
        return [i for i, r in enumerate(ranges) if not (IntRangeSet(r) & irs).isempty]


if __name__ == '__main__':
    import dataset as prd
    from pybedtools import BedTool
    d = prd.Dataset('/groups/price/yakir/data/datasets/1000G3.wim5u/1000G3.wim5u.')
    ldbreakpoints = BedTool(
            '/groups/price/yakir/data/reference/pickrell_breakpoints.hg19.eur.bed')
    sp = SnpPartition(d.ucscbed(22), ldbreakpoints)
    print(sp.ranges())
    print(sp.ind_ranges_overlapping(IntRangeSet('1190:2000')))
