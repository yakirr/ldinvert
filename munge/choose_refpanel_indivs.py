from __future__ import print_function, division
import pyutils.fs as fs
from primitives import *

d = Dataset('GERAhg19')

for nk in [2,3,4,8,10]:
    ai = d.random_indivs_df(nk * 1000)
    print(len(ai))
    dir_name = '/groups/price/yakir/data/datasets/' + d.name + '_ref' + str(nk) + 'k'
    fs.makedir(dir_name)
    ai.to_csv(dir_name + '/indivs_fam', index=False, header=False, sep='\t')
