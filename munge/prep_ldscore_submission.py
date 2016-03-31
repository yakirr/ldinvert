from __future__ import print_function, division
import pandas as pd
from primitives import *

folder = '/groups/price/ldsc/reference_files/UK10K/totalld_allsnps/'

print_filenames = []
extract_filenames = []
chr_filenames = []

with open(folder + 'torun', 'w') as f:
    for chrnum in range(1,23):
        print(chrnum)
        d = Dataset('UK10Khg19.' + str(chrnum))
        print('reading bim')
        d_bim = pd.read_csv(d.genotypes_bedfile.filename + '.bim', sep='\t',
                names=['CHR', 'rsid', 'cM', 'BP', 'A1', 'A2'])
        for i, s in enumerate(d.slices(slice_size=25000)):
            printsnps_filename = '{}.{}.printsnps'.format(d.name, str(i))
            extractsnps_filename = '{}.{}.extractsnps'.format(d.name, str(i))

            b = d.buffer_around_slice(s, 1)
            print(s, b)
            d_bim.iloc[s[0]:s[1], 1].to_csv(
                    folder+printsnps_filename, index=False, header=False)
            d_bim.iloc[b[0]:b[1], 1].to_csv(
                    folder+extractsnps_filename, index=False, header=False)
            print('{}\t{}\t{}'.format(
                printsnps_filename,
                extractsnps_filename,
                chrnum), file=f)
