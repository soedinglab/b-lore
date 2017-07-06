#!/usr/bin/env python

import numpy as np
import os
import ctypes


def czcompgrad(pi, mu, sig2, z, v, nhyp, feat, mureg, reg2, prec, get_gradient=True, is_covariate=False):
    ''' An efficient C code for calculating marginal log likelihood and its gradient if requested.
        This is the only interface with margloglik.so
        The marginal likelihoods are different for: a) regular SNP locus, and b) covariate locus
        Inputs:
              pi: array of N floats
              mu: array of N floats
              sig2: variance of the causal distribution, array of N floats, different for SNP loci and covariate locus
              z: all the zstates for the locus // could be empty for the 0 z-state // all 1's for covariate locus
              v: effect size / coefficients after logistic regression, array of N floats
              nhyp: number of hyperparameters which are being optimized
              feat: the matrix of functional annotations, K x N array
              mureg: the mean of the regularizer, float value
              reg2: the variance of the regularizer, same for SNP loci and covariate locus, float value
              prec: the precision matrix, obtained from the logistic regression, N x N array
              get_gradient: boolean, if false gradient calculation is skipped. Default: True
              is_covariate: boolean, if true calculations are done for a covariate locus
         Returns:
              zcomps: P(z)*F(z) for each zstate
              grad: gradient of each hyperparameter. Returns a zeros array if get_gradient is false
    '''

    _path = os.path.dirname(__file__)
    lib = np.ctypeslib.load_library('../lib/margloglik.so', _path)
    ccomps = lib.margloglik_zcomps
    ccomps.restype = ctypes.c_int
    ccomps.argtypes = [ctypes.c_int,                                                                   # nvar, number of SNPs
                       ctypes.c_int,                                                                   # nhyp, number of hyperparameters
                       ctypes.c_int,                                                                   # nfeat, number of features
                       ctypes.c_int,                                                                   # zlen
                       ctypes.c_int,                                                                   # max(znorm)
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'), # pi
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'), # mu
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'), # sig2
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'), # v^tilde
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'), # Lambda^tilde
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'), # Feature matrix
                       np.ctypeslib.ndpointer(ctypes.c_int, flags='C_CONTIGUOUS, ALIGNED'),            # zarr
                       np.ctypeslib.ndpointer(ctypes.c_int, flags='C_CONTIGUOUS, ALIGNED'),            # znorm
                       ctypes.c_double,                                                                # sigreg2
                       ctypes.c_double,                                                                # mureg
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'), # zcomps
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'), # grad
                       ctypes.c_bool,                                                                  # boolean to get the gradient
                       ctypes.c_bool                                                                   # boolean for the covariate
                      ]

    zlen = len(z)
    zarr = np.array([item for sublist in z for item in sublist], dtype=np.int32)
    znorm = np.array([len(sublist) for sublist in z], dtype=np.int32)

    zcomps = np.zeros(zlen)
    grad = np.zeros(nhyp)

    success = ccomps(v.shape[0],
                     nhyp,
                     feat.shape[0],
                     zlen,
                     np.max(znorm),
                     pi,
                     mu,
                     sig2,
                     v,
                     prec.reshape(-1),
                     feat.reshape(-1),
                     zarr,
                     znorm,
                     reg2,
                     mureg, # mureg
                     zcomps,
                     grad,
                     get_gradient,
                     is_covariate)

    return zcomps, grad


def margloglik_zcomps(pi, mu, sig2, z, v, mureg, reg2, prec, is_covariate=False):
    ''' Calculates the log marginal likelihood for all z-states in a given locus
        Inputs:
              pi: proportion of causal distribution, array of N floats, only applicable for SNP loci
              mu: mean of the causal distribution, array of N floats, 0 for covariate locus
              sig2: variance of the causal distribution, array of N floats, different for SNP loci and covariate locus
              z: all the zstates for the locus // could be empty for the 0 z-state // all 1's for covariate locus
              v: effect size / coefficients after logistic regression, array of N floats
              reg2: the variance of the regularizer, same for SNP loci and covariate locus, float value
              prec: the precision matrix, obtained from the logistic regression, N x N array
              is_covariate: boolean to mark the covariate locus
         Returns:
              zcomps: P(z)*F(z) for each zstate

         Earlier this function used to return 2 other values, which have been depricated.
              !muz: means of the multivariate Gaussian obtained from integration
              !lamz_inv : Inverse of precision matrix of the multivariate Gaussian obtained from integration
    '''

    get_gradient = False
    nhyp = 0
    feat = np.array([])

    zcomps, grad = czcompgrad(pi, mu, sig2, z, v, nhyp, feat, mureg, reg2, prec, get_gradient, is_covariate)

    return zcomps


def mll_grad(pi, mu, sig2, z, v, nhyp, feat, mureg, reg2, prec, is_covariate=False):
    ''' 
        A wrapper for calculating log marginal likelihood and its gradients.
        Includes poor handling of overflow error for very low log marginal likelihoods.
        See czcompgrad for details.
    '''
    get_gradient = True

    zcomps, grad = czcompgrad(pi, mu, sig2, z, v, nhyp, feat, mureg, reg2, prec, get_gradient, is_covariate)
    pl = np.sum(zcomps)
    # To avoid numerical overflow error, negligible marginal likelihoods are rounded off to low values
    if pl < 1e-200:
        print("Warning: Overflow in calculating logarithm of marginal likelihood. Applying constraints")
        pl = 1e-200
    mll = np.log(pl)
    return mll, grad


# ====================================================================================================================================
# Old python implementations. Kept for legacy reasons
#
#def grad_locus_comps(pi, mu, sig2, nhyp, feat, z, zcomps, intmu, intlam_inv, is_covariate=False):
#    ''' Computes the gradient of the log marginal likelihood for each locus.
#        Inputs:
#              pi: proportion of causal distribution, array of N floats, only applicable for SNP loci
#              mu: mean of the causal distribution, array of N floats, 0 for covariate locus
#              sig2: variance of the causal distribution, array of N floats, different for SNP loci and covariate locus
#              feat: array of features, shape K x I for SNP loci, Kc x J for covariate loci
#              nhyp: number of hyperparameters
#              k: number of SNP features, K
#              kc: number of covariate features, Kc
#              z: all the zstates for the locus // could be empty for the 0 z-state // all 1's for covariate locus
#              zcomps: pre-calculated p(z)F(z) for each z-state of the locus
#              intmu: means of the multivariate Gaussian obtained from integration // list of arrays for all zstates
#              intlam_inv: Inverse of precision matrix of the multivariate Gaussian obtained from integration // list of arrays for all zstates
#        Returns:
#              der: The derivates with respect to all hyperparameters, array of (3 x K + Kc) floats
#    '''
#
#    zlen  = len(z)         # Number of zstates for the locus
#    znorm  = [len(z[i]) for i in range(zlen)] # norm of all z-states
#    zmask  = [np.array(z[i]) for i in range(zlen)] # convert to array for masking
#
#    der = np.zeros(nhyp)
#
#    zcomps = np.array(zcomps)
#    pz = zcomps / np.sum(zcomps)
#
#    if not is_covariate: # Calculate values for SNP loci / not for hypothetical covariate locus
#
#        k = feat.shape[0]
#
#        # Partial derivative w.r.t. PI
#        pi_grad = np.zeros((zlen, k))
#        for i in range(zlen):
#            mask = np.zeros(pi.shape)
#            if not znorm[i] == 0:
#                mask[zmask[i]] = 1
#            pi_grad[i, :] = np.einsum('i, ki -> k', mask - pi, feat)
#        der[0:k] = np.einsum('i, ik -> k', pz, pi_grad)
#    
#        # Partial derivative w.r.t. MU
#        mu_grad = np.array([np.einsum('i, ki -> k', \
#                                      (intmu[i] - mu[zmask[i]]) / sig2[zmask[i]], \
#                                      feat[:, zmask[i]]) \
#                           if not znorm[i]==0 else np.zeros(k) for i in range(zlen)])
#        der[k:2*k] = np.einsum('i, ik -> k', pz, mu_grad)
#    
#        # Partial derivative w.r.t. SIGMA
#        a = [(np.square(mu[zmask[i]] - intmu[i]) - sig2[zmask[i]] + np.diag(intlam_inv[i])) / sig2[zmask[i]] \
#                           if not znorm[i]==0 else np.array([]) for i in range(zlen)]    
#        sig_grad = np.array([np.einsum('i, ki -> k', a[i], feat[:, zmask[i]]) \
#                           if not znorm[i]==0 else np.zeros(k) for i in range(zlen)])
#        der[2*k:3*k] = 0.5 * np.einsum('i, ik -> k', pz, sig_grad)
#
#    else:
#
#        kc = feat.shape[0]
#        # Partial derivative w.r.t. SIGMA_W
#        a = [(np.square(mu[zmask[i]] - intmu[i]) - sig2[zmask[i]] + np.diag(intlam_inv[i])) / sig2[zmask[i]] \
#                           if not znorm[i]==0 else np.array([]) for i in range(zlen)]
#        sig_grad = np.array([np.einsum('i, ki -> k', a[i], feat[:, zmask[i]]) \
#                           if not znorm[i]==0 else np.zeros(kc) for i in range(zlen)])
#        der[-kc:] = 0.5 * np.einsum('i, ik -> k', pz, sig_grad)
#
#    return der
#
#
#def margloglik_zcomps_py(pi, mu, sig2, z, v, reg2, prec, is_covariate=False):
#    ''' Calculates the log marginal likelihood for all z-states in a given locus
#        Inputs:
#              pi: proportion of causal distribution, array of N floats, only applicable for SNP loci
#              mu: mean of the causal distribution, array of N floats, 0 for covariate locus
#              sig2: variance of the causal distribution, array of N floats, different for SNP loci and covariate locus
#              z: all the zstates for the locus // could be empty for the 0 z-state // all 1's for covariate locus
#              v: effect size / coefficients after logistic regression, array of N floats
#              reg2: the variance of the regularizer, same for SNP loci and covariate locus, float value
#              prec: the precision matrix, obtained from the logistic regression, N x N array
#              is_covariate: boolean to mark the covariate locus
#         Returns:
#              muz: means of the multivariate Gaussian obtained from integration
#              lamz_inv : Inverse of precision matrix of the multivariate Gaussian obtained from integration
#              zcomps: 
#    '''
#
#    zlen   = len(z) # number of z-states
#    znorm  = [len(z[i]) for i in range(zlen)] # norm of all z-states
#    zmask  = [np.array(z[i]) for i in range(zlen)] # convert to array for masking
#
#    # Log of probability of z given Xi and Theta => ln p(z | Xi, Theta)
#    # For each locus, this is calculated from the array of pi
#    # For covariate locus, p(z | Xi, Theta) is 1, hence ln p is zero.
#    '''Be careful: For any SNP, if pi is zero, then log(pi) becomes undefined
#                   For any SNP, if pi is 1, then log(1-pi)  becomes undefined
#    '''
#
#    logpz = [0.0 for i in range(zlen)] # ln(1) = 0
#    if not is_covariate: # Calculate values for SNP loci / not for hypothetical covariate locus
#        for i in range(zlen):
#            mask = np.zeros(pi.shape, dtype=bool)
#            if not znorm[i] == 0:
#                mask[zmask[i]] = True
#            logpz[i] = np.sum(np.log(pi[mask])) + np.sum(np.log(1.0 - pi[~mask]))
#
#    # Log of second component => ln F(z | Xi, Theta)
#    # First we calculate 2 pre-requisite terms: Lambda_z and V_z
#
#    # Lambda_z
#    # Precision matrix of the multivariate Gaussian obtained from integration.
#    # For locus with ||z|| = 0, it is an empty array
#    # For covariate locus, input sig2 is different, everything else is same
#    def __getlamz(mask):
#        ''' prec, sig2, reg2 comes from global scope '''
#        x = np.array(prec[mask][:,mask])
#        x[np.diag_indices(len(mask), ndim=2)] += 1 / sig2[mask] - 1 / reg2
#        return x
#
#    lamz     = [__getlamz(zmask[i])    if not znorm[i] == 0 else np.array([]) for i in range(zlen)]
#    lamz_inv = [np.linalg.inv(lamz[i]) if not znorm[i] == 0 else np.array([]) for i in range(zlen)]
#
#    # V_z 
#    # Vector of means of the multivariate Gaussian obtained from integration
#    # For locus with ||z|| = 0, it is an empty array
#    # For covariate locus, input mu, sig2 are different, everything else is same
#    vzinner = np.dot(prec, v) + (mu / sig2) # + (0 / reg2) 
#    muz = [np.dot(lamz_inv[i], vzinner[zmask[i]]) if not znorm[i] == 0 else np.array([]) for i in range(zlen)]
#
#    # Finally, we calculate the 4 terms of ln F
#    t1  = [np.sum(np.log(1 / sig2[zmask[i]]))                    if not znorm[i] == 0 else 1 for i in range(zlen)]
#    t2  = [np.linalg.slogdet(lamz[i])                       if not znorm[i] == 0 else (1, -1) for i in range(zlen)]
#    t3  = [np.sum(np.square(mu[zmask[i]]) / sig2[zmask[i]]) if not znorm[i] == 0 else 0 for i in range(zlen)]
#    t4  = [np.einsum('i, ij, j', muz[i], lamz[i], muz[i])   if not znorm[i] == 0 else 0 for i in range(zlen)]
#
#    logfz = [(t2[i][0], 0.5 * (t1[i] - t2[i][1] - t3[i] + t4[i])) for i in range(zlen)]
#
#    zcomps = [np.exp(logpz[i] + logfz[i][1]) for i in range(zlen)]
# 
#    return muz, lamz_inv, zcomps
