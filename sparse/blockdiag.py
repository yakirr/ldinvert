from __future__ import print_function, division
import numpy as np
import copy
import genome.utils as gutils
import pysnptools.util as psutil
from pysnptools.util import IntRangeSet

class BlockDiagArray(object):
    def __init__(self, ranges_to_arrays):
        self.ranges_to_arrays = {}
        for r in ranges_to_arrays.keys():
            if ranges_to_arrays[r].shape and ranges_to_arrays[r].size > 0:
                self.ranges_to_arrays[r] = ranges_to_arrays[r]
        self.intrangeset = IntRangeSet(self.ranges_to_arrays.keys())

    def __str__(self):
        result = ''
        for r, a in self.ranges_to_arrays.items():
            result += str(r) + '\n' + str(a) + '\n'
        return result

    def copy(self):
        return copy.deepcopy(self)

    def ranges(self):
        return self.ranges_to_arrays.keys()

    def inv(self):
        ranges_to_invarrays = {
                r:np.linalg.inv(self.ranges_to_arrays[r]) for r in self.ranges()
                }
        return BlockDiagArray(ranges_to_invarrays)

    def trace(self):
        return np.sum([
            np.trace(a) for a in self.ranges_to_arrays.values()
            ])

    # merge several disjoint BDM's
    @staticmethod
    def merge(bdms):
        ranges = [r for bdm in bdms for r in bdm.ranges()]
        ranges_to_arrays = {
                r:bdm.ranges_to_arrays[r] for bdm in bdms for r in bdm.ranges()
                }
        return BlockDiagArray(ranges_to_arrays)

    # intersect this block with another set of regions
    def __zero_block_outside_irs(self, r, other_intrangeset):
        my_intrangeset = IntRangeSet(r)
        intersection_intrangeset = my_intrangeset & other_intrangeset
        mask = np.zeros(len(my_intrangeset), dtype=bool)
        for s in intersection_intrangeset.ranges():
            start = my_intrangeset.index(s[0])
            end = start + s[1] - s[0]
            mask[start:end] = True
        self.ranges_to_arrays[r][~mask] = 0
        self.ranges_to_arrays[r].T[~mask] = 0 # done this way for compatibility with 1d arrays

    # create a copy of this BDM restricted to a set of ranges.
    def zero_outside_irs(self, other_intrangeset):
        for r in self.ranges():
            self.__zero_block_outside_irs(
                    r,
                    other_intrangeset)
        return self

    # assumes the other matrix has the exact same set of ranges
    def dot(self, other):
        result_ranges_to_arrays = {
                r:self.ranges_to_arrays[r].dot(other.ranges_to_arrays[r])
                for r in self.ranges()
                }
        if result_ranges_to_arrays.values()[0].shape:
            return BlockDiagArray(
                    result_ranges_to_arrays)
        else:
            return sum(result_ranges_to_arrays.values())

    # adds lambda I to the diagonals of everything
    def add_lambdaI(self, Lambda, renormalize=False):
        normalization = 1 + Lambda if renormalize else 1
        return BlockDiagArray({
                r:(self.ranges_to_arrays[r] + Lambda * np.eye(len(A))) / normalization
                for r, A in self.ranges_to_arrays.items()
            })

    # assumes all the arrays are 1-d
    def plot(self, outfile):
        pass #TODO: implement


    # assumes the both arrays have the exact same set of ranges
    @staticmethod
    def solve(A, b):
        result_ranges_to_arrays = {
                r:np.linalg.solve(A.ranges_to_arrays[r], b.ranges_to_arrays[r])
                for r in A.ranges()
                }
        return BlockDiagArray(result_ranges_to_arrays)

    @staticmethod
    def eye(intrangeset):
        return BlockDiagArray({r:np.eye(r[1] - r[0]) for r in intrangeset.ranges()})



# creates a bda from a large (1d) array by restricting it to the relevant intrangeset
def bda_from_big1darray(bigarray, intrangeset):
    return BlockDiagArray({
        r : bigarray[r[0]:r[1]] for r in intrangeset.ranges()})


if __name__ == '__main__':
    from pybedtools import Interval
    np.set_printoptions(precision=3, linewidth=200, suppress=True)
    M = np.random.rand(100, 100)
    intrangeset = IntRangeSet('1:20,50:60,99:100')
    small_regions_to_indexsets = IntRangeSet('3:5,52:54,56:59')

    # create block diagonal matrix manually
    M_blockdiag = np.zeros(M.shape)
    for r in intrangeset.ranges():
        mask = np.zeros(len(M_blockdiag), dtype=bool)
        mask[r[0]:r[1]] = True
        temp = np.copy(M)
        temp[~mask] = 0
        temp.T[~mask] = 0
        M_blockdiag += temp

    # create equivalent BDM
    M_class = BlockDiagArray(
            {r:psutil.sub_matrix(M, range(r[0], r[1]),
                range(r[0], r[1])) for r in intrangeset.ranges()})

    # test inv function of BlockDiagArray
    M_class_copy = M_class.copy()
    I = M_class.dot(M_class_copy.inv())
    print('M_class * M_class_copy^{-1} =', I)

    # test quadratic form evaluation
    # create the vector to operate on together with its BDM version
    v = np.arange(100, dtype=np.float32)
    v_class = bda_from_big1darray(v, intrangeset)

    print('dense vT.M.v =', v.dot(M_blockdiag.dot(v)))
    print('sparse vT.M.V =', v_class.dot(M_class.dot(v_class)))

    # create dense version of the block diagonal matrix restricted to the smaller set
    M_blockdiag_smaller = np.copy(M_blockdiag)
    mask = np.zeros(len(M), dtype=bool)
    mask[small_regions_to_indexsets] = True
    M_blockdiag_smaller[~mask] = 0
    M_blockdiag_smaller.T[~mask] = 0

    # restrict the BDMs to the smaller set
    M_class.zero_outside_irs(small_regions_to_indexsets)
    v_class.zero_outside_irs(small_regions_to_indexsets)

    print('dense vT.M.v =', v.dot(M_blockdiag_smaller.dot(v)))
    print('sparse vT.M.V =', v_class.dot(M_class.dot(v_class)))
