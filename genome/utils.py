from __future__ import print_function, division
import numpy as np
from pybedtools import BedTool, Interval
from pysnptools.snpreader import Bed
from pysnptools.util import IntRangeSet

# restricts a covariance matrix to a class
def make_block_diagonal(R, window_masks):
    result = np.zeros(R.shape)
    for window_mask in window_masks:
        temp = np.copy(R)
        temp[~window_mask] = 0
        temp.T[~window_mask] = 0
        result = result + temp
    return result

# returns length of chromosome
def chromosome_length_in_bp(chrom_num, genome_build='hg19'):
    genome = BedTool('/groups/price/yakir/data/GO.1Mb_flanks/' + genome_build + '.chrlen.bed')
    return int(genome.filter(lambda r: r.chrom == 'chr' + str(chrom_num))[0].end)

# returns the number of SNPs in a given dataset that are in the chromosome
def chromosome_length_in_SNPs(chrom_num, dataset_name):
    chromosome = BedTool('/groups/price/hilary/data/' + dataset_name + '/all.' + str(chrom_num) + '.ucscbed')
    return len(chromosome)

# retrieves the sample size of a given data set
def sample_size(dataset_name):
    dataset = Bed('/groups/price/hilary/data/' + dataset_name + '/all.22')
    return dataset.iid_count

# retrieves the number of snps in a given data set
def snps_in_dataset(dataset_name):
    return sum([chromosome_length_in_SNPs(chrom_num, dataset_name) for chrom_num in range(1,23)])

# return the regions in a region collection that fall in a given chromosome
def regions_in_chromosome(regions, chrom_num):
    return [r for r in regions if r.chrom == 'chr' + str(chrom_num)]

# returns a mask that can be applied to the SNPs in a chromosome
# and pulls out SNPs in a given pathway
def indexset_wrt_chromosome(chrom_num, regions_to_indexsets):
    result = IntRangeSet()
    result.add([regions_to_indexsets[r] for r in regions_to_indexsets.keys() if r.chrom == 'chr' + str(chrom_num)])
    return result

# total number of SNPs in a pathway
def total_size_snps(indexsets):
    return np.sum([len(iset) for iset in indexsets])
