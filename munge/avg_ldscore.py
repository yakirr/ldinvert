from __future__ import print_function, division
import numpy as np

d = Dataset('GERAhg19')
bad_set = SnpSubset(d, GenomicSubset('50').bedtool).expanded_by(0.02)

total = 0
for i in range(500):
    if i in bad_set.irs:
        continue
    m = np.random.randint(0, d.M)
    s = (max(0, m-200), min(d.M, m+200))
    X = d.get_standardized_genotypes(s)
    X1 = d.get_standardized_genotypes((m, m+1))
    corr = X1.T.dot(X) / d.N - 400./d.N
    print(np.linalg.norm(corr)**2)
    total += np.linalg.norm(corr)**2
    print(total / (i+1))
