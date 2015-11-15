from __future__ import print_function, division
import numpy as np
import hyperparams as hp

hp.load()
with hp.results_file() as results:
    header = results.readline().strip()
print(header)

with hp.results_file() as results:
    numbers = np.loadtxt(results, skiprows=1, delimiter='\t')
    num_stats = int((numbers.shape[1] - 2) / 3)

    last_row = np.zeros(numbers.shape[1])
    last_row[0] = last_row[1] = np.nan
    for i in range(num_stats):
        avg_sq_bias = np.var(numbers[:,2+3*i])
        avg_var = np.mean(numbers[:,2+3*i+1])
        avg_mse = np.mean(numbers[:,2+3*i+2])
        last_row[2+3*i : 2+3*i+3] = avg_sq_bias, avg_var, avg_mse

    numbers = np.concatenate((numbers, [last_row]))

with hp.results_file_processed(mode='w') as results:
    np.savetxt(results, numbers, delimiter='\t', header=header, comments='')
