#!/usr/bin/env python

import numpy as np
import os
import ctypes

'''
Granular error handling of margloglik.so
'''
class MarglikError(Exception): pass
class InvalidPiError(MarglikError): pass
class MatrixNotPSDError(MarglikError): pass
class DecompositionError(MarglikError): pass

MARGLIK_ERROR_MAP = {
    1: InvalidPiError,
    2: MatrixNotPSDError,
    3: DecompositionError,
}


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
    ccomps.restype = ctypes.c_double
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
                       ctypes.c_bool,                                                                  # boolean for the covariate
                       ctypes.POINTER(ctypes.c_double)  # logML (output)
                      ]

    zlen = len(z)
    zarr = np.array([item for sublist in z for item in sublist], dtype=np.int32)
    znorm = np.array([len(sublist) for sublist in z], dtype=np.int32)

    logmL = ctypes.c_double()
    zcomps = np.zeros(zlen)
    grad = np.zeros(nhyp)
    retcode = -1

    retcode = ccomps(v.shape[0],
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
                     is_covariate,
                     ctypes.byref(logmL))

    if retcode != 0:
        exc_class = MARGLIK_ERROR_MAP.get(retcode, MarglikError)
        raise exc_class(f"C++ error code {retcode}")

    return logmL.value, zcomps, grad


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
    nhyp = 0 # in C++, nhyp is only required if gradient is calculated.
    feat = np.array([])

    logmL, zprob, grad = czcompgrad(pi, mu, sig2, z, v, nhyp, feat, mureg, reg2, prec, get_gradient, is_covariate)

    return zprob


def mll_grad(pi, mu, sig2, z, v, nhyp, feat, mureg, reg2, prec, is_covariate=False):
    ''' 
        A wrapper for calculating log marginal likelihood and its gradients.
        Includes poor handling of overflow error for very low log marginal likelihoods.
        See czcompgrad for details.
    '''
    get_gradient = True

    logmL, zcomps, grad = czcompgrad(pi, mu, sig2, z, v, nhyp, feat, mureg, reg2, prec, get_gradient, is_covariate)
    return logmL, grad
