plt.close('all')
chrnum_to_beta = pickle.load(hp.beta_file(3))
beta = estimators.to_block_diag(chrnum_to_beta)
regions = gutils.regions_in_chromosome(beta.regions(), 20)
allbeta = np.concatenate([beta.regions_to_arrays[r][np.flatnonzero(beta.regions_to_arrays[r])] for r in regions])

plt.figure()
avgallbetahat = np.zeros(allbeta.shape)
T = 99
for i in range(1,T+1):
    print(i)
    chrnum_to_alphahat = pickle.load(hp.sumstats_file(3, i))
    alphahat = estimators.to_block_diag(chrnum_to_alphahat)
    betahat = bd.BlockDiagArray.solve(estimators.R, alphahat)
    allbetahat = np.concatenate([betahat.regions_to_arrays[r][np.flatnonzero(beta.regions_to_arrays[r])] for r in regions])
    avgallbetahat += allbetahat
    # plt.scatter(allbeta, allbetahat, s=2)
    plt.plot(allbetahat, 'bo')
# plt.scatter(allbeta, allbeta, c=[1,0,0,1])
plt.plot(allbeta, 'ro')
plt.savefig('../plots/beta_betahat.png')

avgallbetahat /= T
plt.figure()
# plt.scatter(allbeta, avgallbetahat)
# plt.scatter(allbeta, allbeta, c=[1,0,0,1])
plt.plot(avgallbetahat, 'bo')
plt.plot(allbeta, 'ro')
plt.title('T=' + str(T))
plt.savefig('../plots/beta_betahat.avg.png')
