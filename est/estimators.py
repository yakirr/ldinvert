from __future__ import print_function, division
import numpy as np
import pickle
import genome.utils as gutils
import sparse.blockdiag as bd
import hyperparams as hp

pathway_regions_to_indexsets = pickle.load(hp.pathway_file())
R = pickle.load(hp.covariance_around_pathway_file())    # R is a BlockDiagArray
R_reg = R.add_lambdaI(0.05, renormalize=True)
R_ridge = R.add_lambdaI(0.05)
RA = R.copy().restrict(pathway_regions_to_indexsets)
RA_Rinv = RA.dot(R.inv())
RA_Rinv_trace = RA_Rinv.trace()
RA_Rreginv_trace = RA.dot(R_reg.inv()).trace()
RA_Rridgeinv_trace = RA.dot(R_ridge.inv()).trace()

print('pathway size:', gutils.total_size_snps(pathway_regions_to_indexsets.values()))
print('merged size:', gutils.total_size_snps(R.indexsets()))

def to_block_diag(chrnum_to_vectors):
    return bd.bda_from_bigarrays(chrnum_to_vectors, R.regions_to_indexsets)

def truth(chrnum_to_beta):
    beta = to_block_diag(chrnum_to_beta)
    return beta.dot(RA.dot(beta))

# alphahat is a BlockDiagArray
def MLE(chrnum_to_alphahat, indiv_indices, Y, N):
    alphahat = to_block_diag(chrnum_to_alphahat)
    betahat = bd.BlockDiagArray.solve(R, alphahat)
    biased = betahat.dot(RA.dot(betahat))
    return biased - RA_Rinv_trace / N

def MLE_reg(chrnum_to_alphahat, indiv_indices, Y, N):
    alphahat = to_block_diag(chrnum_to_alphahat)
    betahat = bd.BlockDiagArray.solve(R_reg, alphahat)
    biased = betahat.dot(RA.dot(betahat))
    return biased - RA_Rreginv_trace / N

def MLE_ridge(chrnum_to_alphahat, indiv_indices, Y, N):
    alphahat = to_block_diag(chrnum_to_alphahat)
    betahat = bd.BlockDiagArray.solve(R_ridge, alphahat)
    biased = betahat.dot(RA.dot(betahat))
    return biased - RA_Rridgeinv_trace / N

def MLE_fixed_indiv(chrnum_to_alphahat, indiv_indices, Y, N):
    alphahat = to_block_diag(chrnum_to_alphahat)
    betahat = bd.BlockDiagArray.solve(R, alphahat)
    biased = betahat.dot(RA.dot(betahat))
    bias = RA_Rinv_trace / N
    return (biased - bias) / (1 - bias)

ldscores_pathway, ldscores_notpathway = hp.ldscores_files()
ldscores_pathway = pickle.load(ldscores_pathway)
ldscores_notpathway = pickle.load(ldscores_notpathway)
keys = ldscores_pathway.keys()
ldscores_pathway = np.concatenate([ldscores_pathway[c] for c in keys])
ldscores_notpathway = np.concatenate([ldscores_notpathway[c] for c in keys])
ldscores = ldscores_pathway + ldscores_notpathway
total_ldscore = np.sum(ldscores)
avg_ldscore = np.mean(ldscores)
X = np.array([np.ones(len(ldscores_pathway)), ldscores_pathway, ldscores_notpathway]).T
print('total ld score to pathway:', np.sum(ldscores_pathway))
print('total ld score to rest of genome:', np.sum(ldscores_notpathway))
def LDSC(chrnum_to_alphahat, indiv_indices, Y, N):
    alphahat = np.concatenate([chrnum_to_alphahat[k] for k in keys])
    chisq = hp.sumstats_params.N * alphahat ** 2

    tauhat = (np.mean(chisq) - 1) / (hp.sumstats_params.N * avg_ldscore)
    w = np.sqrt(np.array([
        1 / (l * (1 + hp.sumstats_params.N * tauhat * l)**2) for l in ldscores
        ]))
    weighted_X = X * w[:,None]
    hat_matrix = np.linalg.inv(weighted_X.T.dot(weighted_X)).dot(weighted_X.T)

    weighted_chisq = w * chisq
    solution = hat_matrix.dot(weighted_chisq)
    result = solution[1] * gutils.total_size_snps(pathway_regions_to_indexsets.values()) / \
            hp.sumstats_params.N
    # import pdb
    # pdb.set_trace()
    return result

methods = {
        'mle' : MLE,
        'mle_reg' : MLE_reg,
        'mle_ridge' : MLE_ridge,
        # 'mle_fixed' : MLE_fixed_indiv,
        'ldsc' : LDSC
        }
