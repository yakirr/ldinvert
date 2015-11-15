from __future__ import print_function, division
import numpy as np
import pickle
import genome.utils as gutils
import sparse.blockdiag as bd
import hyperparams as hp

def global_init(annotation_name):
    global annot_intrangeset, R, RA
    annot_intrangeset = pickle.load(hp.pathway_file(custom_name=annotation_name))
    R = pickle.load(hp.covariance_around_pathway_file(custom_name=annotation_name))
    # TODO: should RA also be regularized? Maybe it doesn't matter as much
    # since it's not being inverted?
    RA = R.copy().zero_outside_irs(annot_intrangeset)
    print('annotation size:', len(annot_intrangeset))
    print('after adding LD bandwidth:', len(R.intrangeset))

def truth(beta):
    beta_bd = bd.bda_from_big1darray(beta, R.intrangeset)
    return beta_bd.dot(RA.dot(beta_bd))

################
## the MLE-based methods
################
def MLE_init(annotation_name):
    global_init(annotation_name)
    global RA_Rinv_trace, Rinv_RA_Rinv
    Rinv = R.inv()
    RA_Rinv = RA.dot(Rinv)
    RA_Rinv_trace = RA_Rinv.trace()
    Rinv_RA_Rinv = Rinv.dot(RA_Rinv)
def MLE(alphahat, indiv_indices, Y, N):
    alphahat_bd = bd.bda_from_big1darray(alphahat, R.intrangeset)
    betahat = bd.BlockDiagArray.solve(R, alphahat_bd)
    biased = betahat.dot(RA.dot(betahat))
    return biased - RA_Rinv_trace / N

def MLE_reg_init(annotation_name, Lambda):
    global_init(annotation_name)
    global Rreg, RA_Rreginv_trace, Rreginv_RA_Rreginv
    Rreg = R.add_lambdaI(Lambda, renormalize=True)
    Rreginv = Rreg.inv()
    RA_Rreginv = RA.dot(Rreginv)
    RA_Rreginv_trace = RA_Rreginv.trace()
    Rreginv_RA_Rreginv = Rreginv.dot(RA_Rreginv)
def MLE_reg(alphahat, indiv_indices, Y, N, Lambda):
    alphahat_bd = bd.bda_from_big1darray(alphahat, Rreg.intrangeset)
    biased = alphahat_bd.dot(Rreginv_RA_Rreginv.dot(alphahat_bd))
    return biased - RA_Rreginv_trace / N

def MLE_ridge_init(annotation_name, Lambda):
    global_init(annotation_name)
    global Rridge, RA_Rridgeinv_trace, Rridgeinv_RA_Rridgeinv
    Rridge = R.add_lambdaI(Lambda)
    Rridgeinv = Rridge.inv()
    RA_Rridgeinv = RA.dot(Rridgeinv)
    RA_Rridgeinv_trace = RA_Rridgeinv.trace()
    Rridgeinv_RA_Rridgeinv = Rridgeinv.dot(RA_Rridgeinv)
def MLE_ridge(alphahat, indiv_indices, Y, N, Lambda):
    alphahat_bd = bd.bda_from_big1darray(alphahat, Rridge.intrangeset)
    biased = alphahat_bd.dot(Rridgeinv_RA_Rridgeinv.dot(alphahat_bd))
    return biased - RA_Rridgeinv_trace / N

def LDSC_init(annotation_name):
    global_init(annotation_name)
    global avg_ldscore, ldscores, X
    print('loading ld scores')
    ldscores_annot = pickle.load(hp.ldscores_file_pathway(custom_name=annotation_name))
    ldscores_notannot = pickle.load(
            hp.ldscores_file_notpathway(custom_name=annotation_name))
    ldscores = ldscores_annot + ldscores_notannot
    total_ldscore = np.sum(ldscores)
    avg_ldscore = np.mean(ldscores)
    X = np.array([np.ones(len(ldscores_annot)), ldscores_annot, ldscores_notannot]).T
    print('total ld score to pathway:', np.sum(ldscores_annot))
    print('total ld score to rest of genome:', np.sum(ldscores_notannot))
def LDSC(alphahat, indiv_indices, Y, N):
    chisq = hp.sumstats_params.N * alphahat ** 2

    tauhat = (np.mean(chisq) - 1) / (hp.sumstats_params.N * avg_ldscore)
    w = np.sqrt(np.array([
        1 / (l * (1 + hp.sumstats_params.N * tauhat * l)**2) for l in ldscores
        ]))
    weighted_X = X * w[:,None]
    hat_matrix = np.linalg.solve(weighted_X.T.dot(weighted_X), weighted_X.T)

    weighted_chisq = w * chisq
    solution = hat_matrix.dot(weighted_chisq)
    result = solution[1] * len(annot_intrangeset) / hp.sumstats_params.N
    return result


def method_with_param(f, Lambda):
    return (lambda alpha_hat, indiv_indices, Y, N :
        f(alpha_hat, indiv_indices, Y, N, Lambda))
def init_with_param(f, Lambda):
    return (lambda annotation_name :
        f(annotation_name, Lambda))
methods = {
    'mle' : method_with_param(MLE_ridge, 0),
    # 'mle_reg,lambda=0.15' : method_with_param(MLE_reg, 0.15),
    # 'mle_ridge,lambda=0.15' : method_with_param(MLE_ridge, 0.15),
    'mle_reg,lambda=0.05' : method_with_param(MLE_reg, 0.05),
    'mle_ridge,lambda=0.05' : method_with_param(MLE_ridge, 0.05),
    'mle_reg,lambda=0.025' : method_with_param(MLE_reg, 0.05),
    'mle_ridge,lambda=0.025' : method_with_param(MLE_ridge, 0.05),
    'ldsc' : LDSC
    }
initializers = {
    'mle' : init_with_param(MLE_ridge_init, 0),
    # 'mle_reg,lambda=0.15' : init_with_param(MLE_reg_init, 0.15),
    # 'mle_ridge,lambda=0.15' : init_with_param(MLE_ridge_init, 0.15),
    'mle_reg,lambda=0.05' : init_with_param(MLE_reg_init, 0.05),
    'mle_ridge,lambda=0.05' : init_with_param(MLE_ridge_init, 0.05),
    'mle_reg,lambda=0.025' : init_with_param(MLE_reg_init, 0.15),
    'mle_ridge,lambda=0.025' : init_with_param(MLE_ridge_init, 0.15),
    'ldsc' : LDSC_init
    }
