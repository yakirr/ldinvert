from __future__ import print_function, division
import math
import numpy as np
from pybedtools import BedTool, Interval
from pysnptools.snpreader import Bed
from pysnptools.util import IntRangeSet
import hyperparams as hp

# returns length of chromosome
def chromosome_length_in_bp(chrom_num, genome_build='hg19'):
    genome = BedTool('/groups/price/yakir/data/reference/' + genome_build + '.chrlen.bed')
    return int(genome.filter(lambda r: r.chrom == 'chr' + str(chrom_num))[0].end)

# returns the number of SNPs in a given dataset that are in the chromosome
#TODO: no longer works because genotypes_bedfile returns the whole genome. needs to be redone.
def chromosome_length_in_SNPs(chrom_num, dataset):
    return dataset.genotypes_bedfile(chrom_num).sid_count

# buffer_size denotes the desired buffer to add to either side of the disjoint slices
def slices(dataset, start=0, end=None, buffer_size=0):
    if not end: end = dataset.M()
    num_slices = int(math.ceil((end - start) / dataset.slice_size))
    for i in xrange(num_slices):
        yield get_slice(dataset, i, start=start, end=end, buffer_size=buffer_size)

def get_slice(dataset, slice_index, start=0, end=None, buffer_size=0):
    if not end: end = dataset.M()
    bufferless_start = start + slice_index * dataset.slice_size
    bufferless_end = min(start + (slice_index+1)*dataset.slice_size, end)
    return (max(bufferless_start - buffer_size, start),
            min(bufferless_end + buffer_size, end))

# # return the regions in a region collection that fall in a given chromosome
# def regions_in_chromosome(regions, chrom_num):
#     return [r for r in regions if r.chrom == 'chr' + str(chrom_num)]

# # returns a mask that can be applied to the SNPs in a chromosome
# # and pulls out SNPs in a given pathway
# def indexset_wrt_chromosome(chrom_num, regions_to_indexsets):
#     result = IntRangeSet()
#     result.add([regions_to_indexsets[r] for r in regions_to_indexsets.keys() if r.chrom == 'chr' + str(chrom_num)])
#     return result
