from __future__ import print_function, division
import numpy as np
import pickle
import genome.utils as gutils
import sparse.blockdiag as bd
import hyperparams as hp

pathway_regions_to_indexsets = pickle.load(hp.pathway_file())
R = pickle.load(hp.covariance_around_pathway_file())    # R is a BlockDiagArray
RA = R.copy().restrict(pathway_regions_to_indexsets)
RA_Rinv = RA.dot(R.inv())
RA_Rinv_trace = RA_Rinv.trace()

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

def MLE_fixed_indiv(chrnum_to_alphahat, indiv_indices, Y, N):
    alphahat = to_block_diag(chrnum_to_alphahat)
    betahat = bd.BlockDiagArray.solve(R, alphahat)
    biased = betahat.dot(RA.dot(betahat))
    bias = RA_Rinv_trace / N
    return (biased - bias) / (1 - bias)

methods = {
        'mle' : MLE,
        'mle_fixed' : MLE_fixed_indiv
        }
