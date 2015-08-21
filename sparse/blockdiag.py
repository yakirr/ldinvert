from __future__ import print_function, division
import numpy as np
import copy
import genome.utils as gutils
import pysnptools.util as psutil
from pysnptools.util import IntRangeSet

class BlockDiagArray(object):
    def __init__(self, regions_to_arrays, regions_to_indexsets):
        self.regions_to_arrays = {}
        self.regions_to_indexsets = {}
        for r in regions_to_arrays.keys():
            if regions_to_arrays[r].shape and regions_to_arrays[r].size > 0:
                self.regions_to_arrays[r] = regions_to_arrays[r]
                self.regions_to_indexsets[r] = regions_to_indexsets[r]

    def __str__(self):
        result = ''
        for r, a in self.regions_to_arrays.items():
            result += str(r) + str(a) + '\n'
        return result

    def copy(self):
        return copy.deepcopy(self)

    def regions(self):
        return self.regions_to_arrays.keys()

    def indexsets(self):
        return self.regions_to_indexsets.values()

    def inv(self):
        regions_to_invarrays = {r:np.linalg.inv(self.regions_to_arrays[r]) for r in self.regions()}
        return BlockDiagArray(regions_to_invarrays, self.regions_to_indexsets)

    def trace(self):
        return np.sum([
            np.trace(a) for a in self.regions_to_arrays.values()
            ])

    # merge several disjoint BDM's
    @staticmethod
    def merge(bdms):
        regions = [r for bdm in bdms for r in bdm.regions()]
        regions_to_arrays = {r:bdm.regions_to_arrays[r] for bdm in bdms for r in bdm.regions()}
        regions_to_indexsets = {r:bdm.regions_to_indexsets[r] for bdm in bdms for r in bdm.regions()}
        return BlockDiagArray(regions_to_arrays, regions_to_indexsets)

    # intersect this block with another set of regions
    @staticmethod
    def __restrict_block(array, indexset, other_indexset):
        intersection_indexset = indexset & other_indexset
        mask = np.zeros(len(indexset), dtype=bool)
        for stretch in intersection_indexset.ranges():
            start = indexset.index(stretch[0])
            mask[start:start + stretch[1] - stretch[0]] = True
        array[~mask] = 0
        array.T[~mask] = 0     # if array is a 1d array then this line will not create a problem

    def restricted_to_chromosome(self, chrom_num):
        regions = gutils.regions_in_chromosome(self.regions(), chrom_num)
        return BlockDiagArray({r:self.regions_to_arrays[r] for r in regions},
                {r:self.regions_to_indexsets[r] for r in regions})

    # create a copy of this BDM restricted to a set of ranges.
    def restrict(self, regions_to_indexsets):
        bdms = []
        for chrom_num in range(1,23):
            self_on_chrom = self.restricted_to_chromosome(chrom_num)
            regions_in_chrom = gutils.regions_in_chromosome(
                    regions_to_indexsets.keys(),
                    chrom_num
                    )
            union_in_chrom_indexset = IntRangeSet()
            union_in_chrom_indexset.add([regions_to_indexsets[r] for r in regions_in_chrom])
            for r in self_on_chrom.regions():
                BlockDiagArray.__restrict_block(
                        self_on_chrom.regions_to_arrays[r],
                        self_on_chrom.regions_to_indexsets[r],
                        union_in_chrom_indexset)
            bdms.append(self_on_chrom)
        return BlockDiagArray.merge(bdms)

    # assumes the other matrix has the exact same set of ranges
    def dot(self, other):
        result_regions_to_arrays = {
                r:self.regions_to_arrays[r].dot(other.regions_to_arrays[r])
                for r in self.regions()
                }
        if result_regions_to_arrays.values()[0].shape:
            return BlockDiagArray(
                    result_regions_to_arrays,
                    self.regions_to_indexsets
                    )
        else:
            return sum(result_regions_to_arrays.values())

    # adds lambda I to the diagonals of everything
    def add_lambdaI(self, Lambda):
        return BlockDiagArray(
                { r:self.regions_to_arrays[r] + Lambda * np.eye(len(A)) for r, A in self.regions_to_arrays.items() },
                self.regions_to_indexsets)

    # assumes the both arrays have the exact same set of ranges
    @staticmethod
    def solve(A, b):
        result_regions_to_arrays = {
                r:np.linalg.solve(A.regions_to_arrays[r], b.regions_to_arrays[r])
                for r in A.regions()
                }
        return BlockDiagArray(
                result_regions_to_arrays,
                A.regions_to_indexsets
                )

    @staticmethod
    def eye(regions_to_indexsets):
        return BlockDiagArray(
                { r:np.eye(len(indexset)) for r, indexset in regions_to_indexsets.items() },
                regions_to_indexsets
                )



def bda_from_bigarrays(chrom_num_to_bigarrays, regions_to_indexsets):
    regions_to_arrays = {}
    for r, indexset in regions_to_indexsets.items():
        chrom_num = int(r.chrom[3:])
        if chrom_num in chrom_num_to_bigarrays:
            regions_to_arrays[r] = chrom_num_to_bigarrays[chrom_num][regions_to_indexsets[r]]
        else:
            regions_to_arrays[r] = np.zeros(len(regions_to_indexsets[r]))
    return BlockDiagArray(regions_to_arrays, regions_to_indexsets)

# assumes that arrays_npz and masks_npz have the same
# set of keys
#TODO: delete the two functions below?
def bdm_from_npz(arrays_npz, masks_npz):
    regions = arrays_npz.keys()
    return BlockDiagArray(
            regions,
            [arrays_npz[r] for r in regions],
            [np.flatnonzero(masks_npz[r]) for r in regions]
            )

def bdm_from_chrom_npz(arrays_npz):
    keys, regions = [], []
    for chrom_num in range(1,23):
        if str(chrom_num) in arrays_npz.keys():
            keys.append(str(chrom_num))
        elif 'chr' + str(chrom_num) in arrays_npz.keys():
            keys.append('chr' + str(chrom_num))
        regions.append('chr' + str(chrom_num) + '\t*\t*\n')
    return BlockDiagArray(
            regions,
            [arrays_npz[k] for k in keys],
            [IntRangeSet((0, len(arrays_npz[k]))) for k in keys]
            )

if __name__ == '__main__':
    from pybedtools import Interval
    np.set_printoptions(precision=3, linewidth=200, suppress=True)
    M = np.random.rand(100, 100)
    regions_to_indexsets = {
            Interval('chr1',1,19):IntRangeSet('1:20'),
            Interval('chr1',50,59):IntRangeSet('50:60'),
            Interval('chr1',99,99):IntRangeSet('99:100')
            }
    small_regions_to_indexsets = {
            Interval('chr1',3,4):IntRangeSet('3:5'),
            Interval('chr1',52,53):IntRangeSet('52:54'),
            Interval('chr1',56,59):IntRangeSet('56:59')
            }

    # create dense block diagonal matrix
    M_blockdiag = np.zeros((100, 100))
    for iset in regions_to_indexsets.values():
        temp = np.copy(M)
        mask = np.zeros(len(temp), dtype=bool)
        mask[iset] = True
        temp[~mask] = 0
        temp.T[~mask] = 0
        M_blockdiag += temp

    # create equivalent BDM
    M_class = BlockDiagArray(
            {r:psutil.sub_matrix(M, iset, iset) for r, iset in regions_to_indexsets.items()},
            regions_to_indexsets
            )

    M_class_copy = M_class.copy()
    I = M_class.dot(M_class_copy.inv())
    print('M_class * M_class_copy^{-1} =', I)

    # create the vector to operate on together with its BDM version
    v = np.arange(100, dtype=np.float32)
    v_class = BlockDiagArray({r:v[iset] for r, iset in regions_to_indexsets.items()}, regions_to_indexsets)

    print('dense vT.M.v =', v.dot(M_blockdiag.dot(v)))
    print('sparse vT.M.V =', v_class.dot(M_class.dot(v_class)))

    # create dense version of the block diagonal matrix restricted to the smaller set
    M_blockdiag_smaller = np.copy(M_blockdiag)
    all_small_indexsets = IntRangeSet.union(*small_regions_to_indexsets.values())
    mask = np.zeros(len(M), dtype=bool)
    mask[all_small_indexsets] = True
    M_blockdiag_smaller[~mask] = 0
    M_blockdiag_smaller.T[~mask] = 0

    # restrict the BDMs to the smaller set
    M_class.restrict(small_regions_to_indexsets)
    v_class.restrict(small_regions_to_indexsets)

    print('dense vT.M.v =', v.dot(M_blockdiag_smaller.dot(v)))
    print('sparse vT.M.V =', v_class.dot(M_class.dot(v_class)))
