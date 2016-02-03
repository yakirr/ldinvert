from __future__ import print_function, division
import copy
import gzip
import numpy as np
import pandas as pd
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

    def restricted_to_chrom_bedtool(self, chrnum):
        return self.bedtool.filter(
                lambda x : x[0] == 'chr' + str(int(chrnum))).saveas()

    @classmethod
    def reference_genome(cls, assembly='hg18'):
        return GenomicSubset(assembly, path=paths.reference, assembly=assembly)

    @classmethod
    def reference_chrom_bedtool(cls, chrnum, assembly='hg18'):
        return cls.reference_genome(assembly=assembly).restricted_to_chrom_bedtool(chrnum)

    @classmethod
    def whole_genome(cls, assembly='hg18'):
        return cls(assembly, path=paths.reference)

class SnpSubset(object):
    def __init__(self, dataset, bedtool=None, irs=None):
        # use bedtools to create an indicator vector for the snps membership in the subset
        self.dataset = dataset
        if bedtool:
            indicator = dataset.snp_coords().intersect(bedtool, c=True)
            self.irs = IntRangeSet(np.flatnonzero(
                np.array([int(snp.name) for snp in indicator])))
        elif irs:
            self.irs = irs
        else:
            self.irs = IntRangeSet()

    def num_snps(self):
        return len(self.irs)

    def expand_by(self, expansion_in_each_direction, units='Morgans'):
        result = IntRangeSet()
        for r in self.irs.ranges():
                result += self.dataset.buffer_around_slice(
                        r, expansion_in_each_direction, units=units)
        self.irs = result

    def expanded_by(self, expansion_in_each_direction, units='Morgans'):
        result = copy.copy(self)
        result.expand_by(expansion_in_each_direction, units=units)
        return result

    # prints subsets in the appropriate format for ldsc
    # all subsets must have the same dataset
    @classmethod
    def print_subsets(cls, outfilename, snpsubsets, names, add_other=False):
        def snp_info_df(d):
            bfile = d.genotypes_bedfile().filename
            return pd.read_csv(bfile + '.bim',
                    delim_whitespace=True,
                    usecols=[0,1,2,3],
                    names=['CHR','SNP','CM','BP'])

        # check that all snpsubsets have the same data set
        if len(set([ss.dataset for ss in snpsubsets])) > 1:
            print('error: all subsets must have the same underlying dataset')
            return
        if not outfilename.endswith('.gz'):
            print('outfilename must end with ".gz". I only write zipped files')
            return

        # get snp info for this dataset
        d = snpsubsets[0].dataset
        df = snp_info_df(d)

        # add the 'other' annotation if necessary
        if add_other:
            union = IntRangeSet()
            for ss in snpsubsets:
                union.update(ss.irs)
            snpsubsets.append(SnpSubset(d, irs=d.all_snps() - union))
            names.append('OTHER')

        # create the pandas dataframe and output it
        for name, ss in zip(names, snpsubsets):
            df[name] = 0
            df.ix[ss.irs, name] = 1
        df = df[['CHR','BP','SNP','CM'] + names]
        with gzip.open(outfilename, 'wt') as write_file:
            df.to_csv(write_file, index=False, sep='\t')


if __name__ == '__main__':
    from dataset import Dataset
    d = Dataset('tinyGERA')
    gs = GenomicSubset('tiny')
    ss = SnpSubset(d, bedtool=gs.bedtool)
    print(ss.num_snps())
    print(GenomicSubset.whole_genome().bedtool)
