from __future__ import print_function, division
import numpy as np
import copy
import math
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
    def ld_matrix(cls, d, snpset_irs, bandwidth, band_units='Morgans', indivs=None,
            output=None, make_consistent_with=None):
        if snpset_irs is None:
            snpset_irs = IntRangeSet((0, d.M))
        if make_consistent_with is None:
            positions_to_flip = np.array([])
        else:
            positions_to_flip = d.snp_consistency_vector(make_consistent_with)
        ranges_to_arrays = {}
        def add_covariance_for_range(r):
            print(r)
            range_size = r[1] - r[0]
            cov = np.zeros((range_size, range_size))
            range_genotypes = d.get_standardized_genotypes(r, indivs=indivs)

            def compute_cov_for_snp(m):
                end = d.buffer_around_snp(m, bandwidth, start=r[0], end=r[1],
                        units=band_units)[1]

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

            # make coding of snps consistent with other dataset
            flip = np.array(IntRangeSet(positions_to_flip) & IntRangeSet((r[0],r[1])),
                    dtype=int) - r[0] # dtype required so we can use empty array as index
            ranges_to_arrays[r][flip] *= -1
            ranges_to_arrays[r][:,flip] *= -1

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
        return sorted([r for r in self.irs.ranges()])

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

    def shape(self):
        return [self.ranges_to_arrays[r].shape for r in self.ranges()]

    def trace(self):
        traces = [np.trace(a) for a in self.ranges_to_arrays.values()]
        for i, r in enumerate(self.ranges_to_arrays.keys()):
            print(r, traces[i], end=', ')
        return np.sum(traces)

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
            for r, A in result_ranges_to_arrays.items():
                print(r, A, end=', ')
            return sum(result_ranges_to_arrays.values())

    # returns a copy in which lambda*I has beed added to everything
    def add_ridge(self, Lambda, renormalize=False):
        normalization = 1 + Lambda if renormalize else 1
        return BlockDiag({
                r:(A + Lambda * np.eye(len(A))) / normalization
                for r, A in self.ranges_to_arrays.items()
            })

    def adjusted_before_inversion(self, Nadjust):
        return BlockDiag({
                r: A * Nadjust / (Nadjust - A.shape[0] - 1.)
                for r, A in self.ranges_to_arrays.items()
            })

    # assumes all the arrays are 2-d. TODO: change that
    def plot(self, irs_to_mark, filename=None):
        import matplotlib.pyplot as plt
        rows = int(math.ceil(math.sqrt(len(self.ranges()))))
        cols = max(int(math.ceil(len(self.ranges()) / rows)), 2) # the max is so we get
                                                                # back a 2-D array always
        fig, axes = plt.subplots(nrows=rows, ncols=cols)
        for i,(r, A) in enumerate(self.ranges_to_arrays.items()):
            width = r[1] - r[0]
            import pdb; pdb.set_trace()
            ax = axes[int(i/rows), i % rows]
            ax.matshow(A, vmin=-1, vmax=1)
            # import pdb; pdb.set_trace()
            my_intrangeset = IntRangeSet(r)
            intersection = my_intrangeset & irs_to_mark
            def draw_line(xs, ys):
                ax.plot(xs, ys, transform=ax.transAxes, lw=0.2, color='k')

            for s in intersection.ranges():
                draw_line([(s[0] - r[0])/width, (s[0] - r[0])/width],
                        [0, 1])
                draw_line([(s[1] - r[0])/width, (s[1] - r[0])/width],
                        [0, 1])
                draw_line([0,1],
                        [(r[1] - s[0])/width, (r[1] - s[0])/width])
                draw_line([0,1],
                        [(r[1] - s[1])/width, (r[1] - s[1])/width])
            ax.set_xticks([0,width]); ax.set_yticks([0,width])
            ax.set_xlim(0,width); ax.set_ylim(width, 0)
            ax.set_title(str(r))

        fig.set_size_inches(axes.shape[0] * 3, axes.shape[1]*4)
        if filename:
            fig.savefig(filename, dpi=400)
        else:
            fig.show()

    # if Nadjust_after is not None, then this bias-adjusts the inverse as if
    # the matrix being inverted is a covariance matrix estimated from a sample of size Nad..
    def inv(self, Nadjust_after=None):
        ranges_to_invarrays = {
                r:np.linalg.inv(self.ranges_to_arrays[r]) for r in self.ranges()
                }
        if Nadjust_after:
            for r in ranges_to_invarrays.keys():
                bias_adjustment = \
                    (Nadjust_after-ranges_to_invarrays[r].shape[0]-1)/Nadjust_after
                ranges_to_invarrays[r] *= bias_adjustment
        return BlockDiag(ranges_to_invarrays)

    # assumes the both arrays have the exact same set of ranges
    # if Nadjust_after is not None, then this bias-adjusts the inverse as if
    # the matrix being inverted is a covariance matrix estimated from a sample of size Nad..
    @staticmethod
    def solve(A, b, Nadjust_after=None):
        result_ranges_to_arrays = {
                r:np.linalg.solve(A.ranges_to_arrays[r], b.ranges_to_arrays[r])
                for r in A.ranges()
                }
        if Nadjust_after:
            for r in result_ranges_to_arrays.keys():
                bias_adjustment = \
                    (Nadjust_after-result_ranges_to_arrays[r].shape[0]-1)/Nadjust_after
                result_ranges_to_arrays[r] *= bias_adjustment
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
