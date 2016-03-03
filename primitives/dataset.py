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
        self.forced_M = forced_M
        self.__dict__.update(
                json.load(open(paths.metadata + 'datasets.json'))[name])
        if chrnum:
            self.name += '|chr' + str(int(chrnum))
            self.bfile += '.' + str(int(chrnum))

        self.__genotypes_bedfile = self.__N = self.__M = self.__chrom_boundaries = None
        self.__covariates = self.__proj = None
        self.__proj_covariates = None

    @property
    def genotypes_bedfile(self):
        if self.__genotypes_bedfile is None:
            if self.forced_M is None:
                self.__genotypes_bedfile = Bed(self.path + self.bfile)
            else:
                self.__genotypes_bedfile = Bed(self.path + self.bfile)[:,:self.forced_M]
        return self.__genotypes_bedfile

    @property
    def N(self):
        if self.__N is None:
            self.__N = self.genotypes_bedfile.iid_count
        return self.__N

    @property
    def M(self):
        if self.__M is None:
            self.__M = self.genotypes_bedfile.sid_count
        return self.__M

    @property
    def chrom_boundaries(self):
        if self.__chrom_boundaries is None:
            # chrom 1 starts at chrom_boundaries[0] and ends before chrom_boundaries[1])
            self.__chrom_boundaries = np.searchsorted(self.genotypes_bedfile.pos[:,0],
                    np.arange(23) +1 )
        return self.__chrom_boundaries

    @property
    def covariates(self):
        if self.__covariates is None:
            print('loading covariates from', self.covariates_file)
            pcs_df = pd.read_csv(self.path + self.covariates_file, sep='\t', header=0)
            fam_df = pd.DataFrame(data=self.genotypes_bedfile.iid, columns=['FID','IID'])
            fam_df['IID'] = fam_df['IID'].astype(int)
            merged_df = pd.merge(fam_df, pcs_df, on='IID').set_index('IID')
            self.__covariates = merged_df.ix[fam_df.ix[:,1], 2:].as_matrix()
        return self.__covariates

    #TODO: this function assumes that once a set of covariates is given it will never change
    #   this can be solved by caching different __proj variables for different covariates
    def project_out_covariates(self, X, covariates=None, constant_term=True):
        if covariates is None:
            new_covariates = self.covariates
        else:
            new_covariates = covariates

        if constant_term:
            new_covariates = np.concatenate(
                (new_covariates, np.ones((new_covariates.shape[0],1))),
                axis=1)

        if self.__proj_covariates is None or \
                not np.allclose(self.__proj_covariates, new_covariates):
            print('computing hat matrix for covariate projection')
            self.__proj_covariates = new_covariates
            self.__proj = np.linalg.solve(
                    self.__proj_covariates.T.dot(self.__proj_covariates),
                    self.__proj_covariates.T)
        result = X - self.__proj_covariates.dot(self.__proj.dot(X))
        print('variance before projection:', np.var(X, axis=0).reshape((-1,))[:10])
        print('variance after projection:', np.var(result, axis=0).reshape((-1,))[:10])
        return result

    def snp_coords(self):
        ucscbedfilename = self.path + self.bfile + '.ucscbed'
        if os.path.exists(ucscbedfilename):
            return BedTool(self.path + self.bfile + '.ucscbed')
        else:
            print('warning: ucscbedfile not found:', ucscbedfilename)
            return None

    def mhc_bedtool(self):
        mhc_filename = paths.reference + self.reference_genome + '.MHC.bed'
        return BedTool(mhc_filename)

    def all_snps(self):
        return IntRangeSet((0, self.M))

    def chromosomes(self):
        return np.unique(self.genotypes_bedfile.pos[:,0]).astype(np.int16)

    def snp_at_distance(self, snp_index, distance, break_ties_to_right=True,
            units='Morgans'):
        if units == 'SNPs':
            result = snp_index + distance
            return max(0, min(result, self.M-1))
        if units == 'Morgans':
            i = 1
        elif units == 'bp':
            i = 2
        else:
            print('ERROR: units must be either Morgans or bp')
            return None

        snp_index = min(snp_index, self.M-1)
        coords = self.genotypes_bedfile.pos[:,i]
        chrom = self.genotypes_bedfile.pos[snp_index,0]
        chrom_coords = coords[self.chrom_boundaries[chrom-1]:self.chrom_boundaries[chrom]]

        if break_ties_to_right:
            return bisect.bisect_right(chrom_coords,
                    coords[snp_index] + distance) - 1 + self.chrom_boundaries[chrom-1]
        else:
            return bisect.bisect_left(chrom_coords,
                    coords[snp_index] + distance) + self.chrom_boundaries[chrom-1]

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
            genotypes = self.genotypes_bedfile[:,s[0]:s[1]].read()
        else:
            genotypes = self.genotypes_bedfile[indivs, s[0]:s[1]].read()
        genotypes.standardize(); genotypes.standardize()
        return genotypes.val

    def random_indivs(self, Nref, replace=False):
        return np.random.choice(self.N, size=Nref, replace=replace)

    def random_indivs_df(self, Nref, replace=False):
        indivs = self.random_indivs(Nref, replace=replace)
        all_indivs = pd.read_csv(self.genotypes_bedfile.filename + '.fam',
                delim_whitespace=True,
                header=None,
                usecols=[0,1])
        return all_indivs.ix[indivs]

    # returns a vector of indices of snps in this dataset that need to have their
    # sign flipped in order to make this dataset consistent with other
    def snp_consistency_vector(self, other):
        my_bimfile = self.genotypes_bedfile.filename + '.bim'
        other_bimfile = other.genotypes_bedfile.filename + '.bim'

        pos_other = pd.read_csv(other_bimfile, sep='\t',
                names=['chr', 'rsid', 'Morgans', 'bp', 'A1', 'A2'])
        pos_me = pd.read_csv(my_bimfile, sep='\t',
                names=['chr', 'rsid', 'Morgans', 'bp', 'A1', 'A2'])
        pos_me.set_index('rsid', append=True)

        merged = pos_me.reset_index().merge(
                pos_other, how='left', on=['rsid']).set_index('index')
        return merged.loc[(merged['A1_x'] == merged['A2_y']) &
                                (merged['A2_x'] == merged['A1_y'])].index.values


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
        print(d.genotypes_bedfile.pos[bs[0],2], d.genotypes_bedfile.pos[s[0],2])
