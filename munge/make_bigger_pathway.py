from __future__ import print_function, division
import os
import numpy as np
import pickle
from pybedtools import BedTool, Interval
from pysnptools.util.intrangeset import IntRangeSet
import pyutils.fs as fs
import hyperparams as hp
import genome.utils as gutils

hp.load()

suffix = '_bigger'
directory = '/groups/price/yakir/data/GO/'
amount_to_add_Kb = 15

# Load up pathway
print('loading pathway')
pathway = BedTool(directory + hp.pathway.name + '.bed').sort()

print(pathway)
print('total size before window addition:', pathway.total_coverage(), 'bp')

# Intersect the pathway with the appropriate genome build
# TODO: this step should be unnecessary if the pathways are correct
print('removing regions not in reference')
chromosomes = BedTool(hp.paths.reference + hp.dataset.reference_genome + '.chrlen.bed')
pathway = chromosomes.intersect(pathway).sort().saveas()

# compute the flanks
print('computing flanks')
flanks = pathway.flank(
    genome=hp.dataset.reference_genome,
    b=amount_to_add_Kb*1000).sort().merge().saveas()
print(flanks)

# compute the union of the flanks and the pathway
print('computing union')
union = pathway.cat(flanks, postmerge=False).sort()
merged = union.merge().saveas(
    directory + hp.pathway.name + suffix + '.bed')
print(merged)
print('total size after window addition:', merged.total_coverage(), 'bp')
