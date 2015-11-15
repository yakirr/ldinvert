from __future__ import print_function, division
import os
import numpy as np
import pickle
from pybedtools import BedTool, Interval
from pysnptools.util.intrangeset import IntRangeSet
import pyutils.fs as fs
import hyperparams as hp
import genome.utils as gutils

def save_snpbelonging(bedfile, outfile):
    snps_ucscbed = hp.dataset.snps_bedtool()
    pathway_indicator = snps_ucscbed.intersect(bedfile, c=True) # create indicator vector
    pathway_indices = IntRangeSet(np.flatnonzero(
            np.array(
                [int(snp.name) for snp in pathway_indicator]
            )))
    print('total size of', outfile.name, '=', len(pathway_indices), 'SNPs')
    pickle.dump(pathway_indices, outfile, 2)

hp.load()

# Load up pathway
print('loading pathway')
pathway = BedTool('/groups/price/yakir/data/GO/' + hp.pathway.name + '.bed').sort()

print(pathway)
print('total size before window addition:', pathway.total_coverage(), 'bp')

# Create directory for output
fs.makedir(hp.paths.pathway_details)

# Intersect the pathway with the appropriate genome build
# TODO: this step should be unnecessary if the pathways are correct
print('removing regions not in reference')
chromosomes = BedTool(hp.paths.reference + hp.dataset.reference_genome + '.chrlen.bed')
pathway = chromosomes.intersect(pathway).sort().saveas(
    hp.pathway_ucscbedfilename())

# compute the flanks
# TODO: use 1cM instead of 1Mb
print('computing flanks')
flanks = pathway.flank(
        genome=hp.dataset.reference_genome,
        b=hp.pathway.window_size_Mb*1000000).sort().merge().saveas(
    hp.pathway_flanks_ucscbedfilename())
print(flanks)

# compute the union of the flanks and the pathway
print('computing union')
union = pathway.cat(flanks, postmerge=False).sort()
merged = union.merge().saveas(
    hp.pathway_with_flanks_ucscbedfilename())
print(merged)
print('total size after window addition:', merged.total_coverage(), 'bp')

# create an indicator vector for whether each SNP in the relevant
# data set is in the pathway
print('computing SNP indexset for each region')
save_snpbelonging(pathway, hp.pathway_file(mode='wb'))
save_snpbelonging(flanks, hp.pathway_flanks_file(mode='wb'))
save_snpbelonging(merged, hp.pathway_with_flanks_file(mode='wb'))
