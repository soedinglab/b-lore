#!/usr/bin/env python

import numpy as np
import collections
import time

from inference.empirical_bayes import EmpiricalBayes
from inference import quasi_laplace

from iotools.io_summary import ReadSummary
from iotools.io_metagwas import WriteResult
from inference.summary_stat import SummaryStatistics
from snptools.features import FunctionalAnnotations
from utils import hyperparameters
from utils.result import create_locus_result


def optimize_hyperparameters(input_files, feature_files, outdir, file_prefix, cmax, muvar, params):

    '''Main calling function of the optimization.
       It performs the following operations:
           * read and combine the input summary statistics
           * read the functional annotations if any
           * run the empirical Bayes optimization
           * calculate causal probability of each SNP and each locus
           * write the result
    '''

    start_time = time.time()

    # Read and combine the input summary statistics
    summary_stat = SummaryStatistics(input_files)
    summary_stat.combine()
    snpinfo  = summary_stat.snpinfo
    vmin     = summary_stat.vmin
    precll   = summary_stat.precll
    sigreg   = summary_stat.sigreg
    mureg    = summary_stat.mureg
    iscov    = summary_stat.iscov
    locnames = summary_stat.locnames

    nloci = len(snpinfo)
    nsnps = [len(x) for x in snpinfo]
    sigreg2 = sigreg * sigreg

    func_annot = FunctionalAnnotations(snpinfo, iscov, feature_files)
    features = func_annot.features

    preprocend_time = time.time()


    # Run the empirical Bayes optimization
    emp_bayes = EmpiricalBayes(vmin, precll, features, iscov, mureg, sigreg2, 1, muvar, params = params)
    emp_bayes.fit()
    params = emp_bayes.params
    if cmax > 1:
        emp_bayes = EmpiricalBayes(vmin, precll, features, iscov, mureg, sigreg2, cmax, False, params = params, rerun = True)
        emp_bayes.fit()
    params = emp_bayes.params
    zstates = emp_bayes.zstates

    # Get the final results
    # There are two main outputs: a) probability of a locus to be causal
    #                             b) probability of a SNP to be causal
    # -------------------------------------------------------------------
    resultlist = list()

    for l in range(nloci):

        pi, mu, sig2 = hyperparameters.transform(params, features[l], iscov[l])
        zcomps = quasi_laplace.margloglik_zcomps(pi, mu, sig2, zstates[l], vmin[l], mureg, sigreg2, precll[l], iscov[l])
        locusres = create_locus_result(locnames[l], iscov[l], snpinfo[l], vmin[l], precll[l], zstates[l], zcomps, pi, mu, sig2)
        resultlist.append(locusres)

    optimend_time = time.time()

    # Write output
    output = WriteResult(outdir, file_prefix)
    output.set_result(resultlist, sigreg, mureg)
    output.write()

    end_time = time.time()

    # Log the time taken
    total_time = end_time - start_time
    read_time = readend_time - start_time
    preproc_time = preprocend_time - readend_time
    optim_time = optimend_time - preprocend_time
    write_time = end_time - optimend_time
    output.write_time(total_time, preproc_time, optim_time, write_time)
