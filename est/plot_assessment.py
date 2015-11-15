from __future__ import print_function, division
import numpy as np
import hyperparams as hp
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

hp.load()
with hp.results_file() as results:
    header = results.readline().strip().split()
header = header[2:]
print(header)

with hp.results_file() as results:
    numbers = np.loadtxt(results, skiprows=1, delimiter='\t')
    last_row = numbers[-1,2:]

    print(last_row)
    plt.plot(range(len(last_row)), last_row)
    plt.xticks(range(len(last_row)), header)
    plt.show()

