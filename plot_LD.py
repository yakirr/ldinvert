r = gutils.regions_in_chromosome(estimators.R.regions(), 20)[0]
r_start = estimators.R.regions_to_indexsets[r].min()
to_highlight = gutils.regions_in_chromosome(estimators.pathway_regions_to_indexsets, 20)

plt.close('all')
plt.figure()
A = estimators.R.regions_to_arrays[r].copy()
for s in to_highlight:
    for start, stop in estimators.pathway_regions_to_indexsets[s].ranges():
        start -= r_start; stop -= r_start
        A[start-1,:] = A[:, start-1] = -1
        A[stop,:] = A[:, stop] = -1
        print(A[start:stop, start:stop])
plt.colorbar(plt.imshow(A, interpolation=None))


plt.savefig('../plots/ld.png')
