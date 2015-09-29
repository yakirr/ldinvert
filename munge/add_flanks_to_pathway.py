from __future__ import print_function, division
import os
import numpy as np
import pickle
from pybedtools import BedTool, Interval
from pysnptools.util.intrangeset import IntRangeSet
import pyutils.fs as fs
import hyperparams as hp
import genome.utils as gutils

import pdb
def snp_belonging(region):
    snps_in_chrom = hp.dataset.snps_bedtool(region.chrom[3:])
    indicator = snps_in_chrom.intersect(BedTool([region]), c=True)
    return IntRangeSet(np.flatnonzero(
            np.array(
                [int(snp.name) for snp in indicator]
            )))

def save_snpbelonging(bedfile, outfile):
    result = {}
    for region in bedfile:
        new_region = Interval(region.chrom, region.start, region.end)
        result[new_region] = snp_belonging(new_region)
    print('total size of', outfile.name, '=', gutils.total_size_snps(result.values()), 'SNPs')
    pickle.dump(result, outfile, 2)

# Load up pathway
print('loading pathway')
pathway = BedTool('/groups/price/yakir/data/GO/' + hp.pathway.name + '.bed').sort()

# These two lines causes the script to restrict the pathway to chromosome 22
if hp.pathway.chr22only:
    pathway = pathway.filter(lambda x: x.chrom == 'chr22').saveas()

print(pathway)
print('total size before window addition:', pathway.total_coverage(), 'bp')

# Create directory for output
fs.makedir(hp.paths.pathway_details)

# Intersect the pathway with the appropriate genome build
# TODO: this step should be unnecessary if the pathways are correct
print('removing regions not in reference')
chromosomes = BedTool(hp.paths.reference + hp.dataset.reference_genome + '.chrlen.bed')
pathway = chromosomes.intersect(pathway).sort().saveas(
        hp.paths.pathway_details + 'pathway.bed')

# compute the flanks
# TODO: use 1cM instead of 1Mb
print('computing flanks')
flanks = pathway.flank(
        genome=hp.dataset.reference_genome,
        b=hp.pathway.window_size_Mb*1000000).sort().merge().saveas(
        hp.paths.pathway_details + 'flanks.bed')
print(flanks)

# compute the union of the flanks and the pathway
print('computing union')
union = pathway.cat(flanks, postmerge=False).sort()
merged = union.merge().saveas(
        hp.paths.pathway_details + 'merged.bed')
print(merged)
print('total size after window addition:', merged.total_coverage(), 'bp')

# create an indicator vector for whether each SNP in the relevant
# data set is in the pathway
print('computing SNP indexset for each region')
save_snpbelonging(pathway, hp.pathway_file(mode='wb'))
save_snpbelonging(flanks, hp.pathway_flanks_file(mode='wb'))
save_snpbelonging(merged, hp.pathway_with_flanks_file(mode='wb'))
