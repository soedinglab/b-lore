#!/usr/bin/env python

import numpy as np
import collections

RESULT_FIELDS = ['locusname', 'iscov', 'causal_prob', 'finemap', 'zstates', 'hessian']
class ResultObj(collections.namedtuple('_ResultObj', RESULT_FIELDS)):
    __slots__ = ()

    @property
    def nsnps(self):
        return len(self.finemap)


SNPRESULT_FIELDS = ['rsid', 'bp_location', 'ref_allele', 'alt_allele', 'vmin', 'prob', 'pi', 'mu', 'sigma']
class SnpResultObj(collections.namedtuple('_SnpResultObj', SNPRESULT_FIELDS)):
    __slots__ = ()


ZSTATERESULT_FIELDS = ['index', 'prob', 'zcomp']
class ZstateResultObj(collections.namedtuple('_ZstateResultObj', ZSTATERESULT_FIELDS)):
    __slots__ = ()


def create_locus_result(name, iscov, snpinfo, vmin, hess, zstates, zcomps, snp_pi, snp_mu, snp_sigma2):

    # zstate results
    zprob = zcomps / np.sum(zcomps)
    zstateres = list()
    for i, z in enumerate(zstates):
        res = ZstateResultObj(index = z, prob = zprob[i], zcomp = zcomps[i])
        zstateres.append(res)

    # finemapping results
    causal_prob = np.zeros(len(snpinfo))
    for i, z in enumerate(zstates):
        for snp in z:
            causal_prob[snp] += zprob[i]

    snpres = list()
    for i, snp in enumerate(snpinfo):
        res = SnpResultObj(rsid = snp.rsid,
                           bp_location = snp.bp_location,
                           ref_allele = snp.ref_allele,
                           alt_allele = snp.alt_allele,
                           vmin = vmin[i],
                           prob = causal_prob[i],
                           pi = snp_pi[i],
                           mu = snp_mu[i],
                           sigma = np.sqrt(snp_sigma2[i]))
        snpres.append(res)

    # Main locus result object
    if iscov:
        locusprob = zprob[0]
    else:
        locusprob = 1 - zprob[0]
    locusres = ResultObj(locusname = name,
                         iscov = iscov,
                         causal_prob = locusprob,
                         finemap = snpres,
                         zstates = zstateres,
                         hessian = hess)

    return locusres
