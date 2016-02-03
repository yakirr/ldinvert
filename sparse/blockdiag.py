from __future__ import print_function, division
import numpy as np
import copy
import pysnptools.util as psutil
from pysnptools.util import IntRangeSet
from pyutils import iter as it

class BlockDiag(object):
    def __init__(self, ranges_to_arrays, irs=None):
        self.ranges_to_arrays = ranges_to_arrays
        if irs is None:
            irs = IntRangeSet(ranges_to_arrays.keys())
        self.irs = irs

    @classmethod
    def from_big1darray(cls, bigarray, irs):
        return cls({r : bigarray[r[0]:r[1]] for r in irs.ranges()}, irs)

    @classmethod
    def ld_matrix(cls, d, snpset_irs, bandwidth_morg, indivs=None, output=None):
        if snpset_irs is None:
            snpset_irs = IntRangeSet((0, d.M))
        ranges_to_arrays = {}
        def add_covariance_for_range(r):
            print(r)
            range_size = r[1] - r[0]
            cov = np.zeros((range_size, range_size))
            range_genotypes = d.get_standardized_genotypes(r, indivs=indivs)

            def compute_cov_for_snp(m):
                end = d.buffer_around_snp(m, bandwidth_morg, start=r[0], end=r[1])[1]

                window_start = m - r[0]
                window_end = end - r[0]
                window = range_genotypes[:, window_start:window_end]

                cov_to_snps_in_window = \
                        range_genotypes[:,m-r[0]].T.dot(window) / range_genotypes.shape[0]
                cov_to_snps_in_window[0] /= 2 # since we're going to symmetrize later

                cov[m-r[0], window_start:window_end] = cov_to_snps_in_window
            map(compute_cov_for_snp, it.show_progress(range(r[0], r[1])))

            # symmetrization
            ranges_to_arrays[r] = cov + cov.T

        map(add_covariance_for_range,
                snpset_irs.ranges())
        return cls(ranges_to_arrays, snpset_irs)

    def __str__(self):
        result = ''
        for r, a in self.ranges_to_arrays.items():
            result += str(r) + '\n' + str(a) + '\n'
        return result

    def copy(self):
        return copy.deepcopy(self)

    def ranges(self):
        return self.irs.ranges()

    def to_dense(self, ignore_ranges=True):
        if not ignore_ranges:
            pass #TODO
        else:
            dense_shape = np.sum(self.shape(), axis=0)
            result = np.zeros(dense_shape)
            curr_dim = 0
            for r in self.ranges():
                length = r[1] - r[0]
                if len(dense_shape) == 2:
                    result[curr_dim:curr_dim+length, curr_dim:curr_dim+length] = \
                            self.ranges_to_arrays[r]
                elif len(dense_shape) == 1:
                    result[curr_dim:curr_dim+length] = self.ranges_to_arrays[r]
                curr_dim += length
        return result

    def inv(self):
        ranges_to_invarrays = {
                r:np.linalg.inv(self.ranges_to_arrays[r]) for r in self.ranges()
                }
        return BlockDiagArray(ranges_to_invarrays)

    def shape(self):
        return [self.ranges_to_arrays[r].shape for r in self.ranges()]

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

    # zero out this BDM outside a set of ranges
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
            return BlockDiag(
                    result_ranges_to_arrays)
        else:
            return sum(result_ranges_to_arrays.values())

    # adds lambda I to the diagonals of everything
    def add_ridge(self, Lambda, renormalize=False):
        normalization = 1 + Lambda if renormalize else 1
        return BlockDiag({
                r:(self.ranges_to_arrays[r] + Lambda * np.eye(len(A))) / normalization
                for r, A in self.ranges_to_arrays.items()
            })

    # assumes all the arrays are 2-d. TODO: change that
    def plot(self, outfile):
        import matplotlib
        pass #TODO: implement


    # assumes the both arrays have the exact same set of ranges
    @staticmethod
    def solve(A, b):
        result_ranges_to_arrays = {
                r:np.linalg.solve(A.ranges_to_arrays[r], b.ranges_to_arrays[r])
                for r in A.ranges()
                }
        return BlockDiag(result_ranges_to_arrays)

    @staticmethod
    def eye(irs):
        return BlockDiag({r:np.eye(r[1] - r[0]) for r in irs.ranges()})


if __name__ == '__main__':
    from primitives import Dataset, GenomicSubset, SnpSubset
    import copy
    from time import time
    import argparse
    np.random.seed(0)
    parser = argparse.ArgumentParser()
    parser.add_argument('--M', type=int, required=True, help='the number of SNPs to use')
    parser.add_argument('-check_dense', action='store_true', default=False)
    args = parser.parse_args()

    d = Dataset('GERA', forced_M=args.M)
    indivs = d.random_indivs(200)
    tiny_gs = GenomicSubset('50')
    tiny_ss = SnpSubset(d, bedtool=tiny_gs.bedtool)
    tiny_buffered_ss = tiny_ss.expanded_by(0.01)

    t0 = time()
    R = BlockDiag.ld_matrix(d, tiny_buffered_ss.irs, 0.01, indivs=indivs) # 1 cM bandwidth
    R = R.add_ridge(0.05, renormalize=True)
    print('trace of renormalized R should be close to M (with noise due to sample vs pop LD',
            R.trace(), tiny_buffered_ss.num_snps(),
            R.trace() == tiny_buffered_ss.num_snps())
    print('computing R took', time() - t0)
    print('shape of R is:', R.shape())

    RA = R.copy()
    RA.zero_outside_irs(tiny_ss.irs)
    b = BlockDiag.from_big1darray(np.random.randn(d.M), R.irs)

    # check inverse computation
    t0 = time()
    Rinvb = BlockDiag.solve(R, b)
    print('R^{-1}b took', time() - t0)
    if args.check_dense:
        Rinvb_dense = np.linalg.solve(R.to_dense(), b.to_dense())
        print('R^{-1}b behaves well:', np.allclose(Rinvb.to_dense(), Rinvb_dense))
        print('distance from dense solution:', np.linalg.norm(Rinvb.to_dense() - Rinvb_dense))

    t0 = time()
    TrRinvRA = BlockDiag.solve(R, RA).trace()
    print('Tr(Rinv*RA) took', time() - t0)
    if args.check_dense:
        TrRinvRA_dense = np.trace(np.linalg.solve(R.to_dense(), RA.to_dense()))
        print('Tr(R*RA) behaves well:', np.allclose(TrRinvRA, TrRinvRA_dense))
        print(TrRinvRA, 'vs', TrRinvRA_dense)
