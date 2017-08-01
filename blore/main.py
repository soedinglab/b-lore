#!/usr/bin/env python

#############################################################################
######                                                                 ######
######         BAYESIAN MULTIVARIATE META-ANALYSIS OF GWAS             ######
######                                                                 ######
#############################################################################
######                                                                 ######
######                           v0.1.2                                ######
######                                                                 ######
######            developed on Ubuntu 15.04 Python 3.4                 ######
######                                                                 ######
#############################################################################
#############################################################################
######                                                                 ######
######  A Bayesian logistic regression approach to                     ######
######      a. create summary statistics from loci data                ######
######      b. finemap SNPs and predict causal loci                    ######
######  in a GWA study.                                                ######
######                                                                 ######
#############################################################################

import numpy as np
import argparse
import os

import individual
import meta

from utils import hyperparameters

def parse_args():
    ''' Parse the arguments passed to the program.
        The program can be run with two options:
            a) --summary : creates summary statistics from genotype and phenotype data
            b) --meta : performs meta-analysis using the multivariate summary statistics
        Options available when using --summary are as follows:
            --gen
            --sample
            --pheno
            --regoptim
            --out
            --prefix
            --reg
            --pca
        Options available when using --meta are as follows:
            --statinfo
            --feature
            --out
            --prefix
            --zmax
            --params
            --muvar
    '''

    parser = argparse.ArgumentParser(description='Bayesian logistic regression for multivariate meta-analysis')

    parser.add_argument('--summary',
                        dest='summary',
                        action='store_true',
                        help='creates summary statistics from individual level phenotype and genotype')

    parser.add_argument('--meta',
                        dest='meta',
                        action='store_true',
                        help='performs meta-analysis using the multivariate summary statistics')
    
    parser.add_argument('--gen',
                        nargs='*',
                        type=str,
                        dest='genotypefiles',
                        metavar='FILE',
                        help='input genotype file(s) -- separate files for separate locus')
    
    parser.add_argument('--sample',
                        type=str,
                        dest='samplefile',
                        metavar='FILE',
                        help='input sample file')

    parser.add_argument('--pheno',
                        default='pheno',
                        type=str,
                        dest='pheno',
                        metavar='STR',
                        help='name of the phenotype as it appears in sample file')

    parser.add_argument('--statinfo',
                        nargs='*',
                        type=str,
                        dest='statfiles',
                        metavar='FILE',
                        help='input file prefixes of summary statistics for different studies to be used for meta-analysis')

    parser.add_argument('--feature',
                        nargs='*',
                        type=str,
                        dest='featurefiles',
                        metavar='FILE',
                        help='input files for functional annotations')
    
    parser.add_argument('--out',
                        dest='outdir',
                        metavar='FILE',
                        help='output directory')

    parser.add_argument('--prefix',
                        dest='file_prefix',
                        metavar='FILE',
                        help='prefix for the output files')

    parser.add_argument('--regoptim',
                        dest='regoptim',
                        action='store_true',
                        help='optimize the regularizer')

    parser.add_argument('--reg',
                        default=0.1,
                        type=float,
                        dest='sigreg',
                        metavar='REAL',
                        help='variance (sigma) for the regularising Gaussian')

    parser.add_argument('--pca',
                        default=0,
                        type=int,
                        dest='npca',
                        metavar='INT',
                        help='Number of Principal Components to be calculated')


    parser.add_argument('--cov',
                        nargs='*',
                        type=str,
                        dest='covariates',
                        metavar='STR',
                        help='names of covariates to use from the sample file')


    parser.add_argument('--zmax',
                        default=2,
                        type=int,
                        dest='zmax',
                        metavar='INT',
                        help='maximum number of iterations for zstates')

    parser.add_argument('--params',
                        nargs='*',
                        type=float,
                        dest='params',
                        metavar='FLOAT',
                        help='initialization parameters [pi, mu, sigma, sigma_cov]')

    parser.add_argument('--muvar',
                        dest='muvar',
                        action='store_true',
                        help='whether there will be any bounds for mu')


    
    opts = parser.parse_args()
    return opts


class InputError(Exception):
    ''' Raise when incorrect options are used in the argument '''
    pass

# ==============================================================================
# Main program
# ==============================================================================

# Get the arguments
opts = parse_args()


# Check if either --summary of --meta is specified
try:
    assert (opts.meta or opts.summary)
except AssertionError:
    print ('Input error: You must specify either --meta or --summary. See --help for details.')
    raise

# while --summary and --meta cannot be specified together
try:
    assert opts.meta != opts.summary
except AssertionError:
    print ('Input error: You cannot specify both --meta and --summary together. See --help for details.')
    raise

if (opts.summary):
    # check for genotypefiles and samplefile if --summary
    try:
        assert (opts.genotypefiles is not None) and (opts.samplefile is not None)
    except AssertionError:
        print('Input error: You must specify --genotype and --phenotype when using with --summary. See --help for details')
        raise

    genotypefiles = opts.genotypefiles
    samplefile = opts.samplefile
    phenotype = opts.pheno
    covnames = opts.covariates
    sigreg = opts.sigreg
    mureg = 0.0
    tolerance = 0.0001
    npca = opts.npca
    regoptim = opts.regoptim

    # check for at least 1 genotypefile
    try:
        assert len(genotypefiles) > 0
    except AssertionError:
        print('Input error: Please specify genotype files')
        raise

    # check that all genotypefiles exist
    for genotypefile in genotypefiles:
        try:
            assert os.path.isfile(genotypefile)
        except AssertionError:
            print('Input error: %s does not exist' % genotypefile)
            raise

    # check that samplefile exists
    try:
        assert os.path.isfile(samplefile)
    except AssertionError:
        print('Input error: Please specify a sample file')
        raise

    if opts.outdir is not None:
        outdir = os.path.realpath(opts.outdir)
    else:
        outdir = os.path.dirname(os.path.realpath(genotypefiles[0]))
    if opts.file_prefix is not None:
        file_prefix = opts.file_prefix
    else:
        studyname = os.path.split(outdir)[1]
        file_prefix = '%s_summary' % studyname

    individual.summary(samplefile, genotypefiles, phenotype, covnames, mureg, sigreg, tolerance, npca, outdir, file_prefix, regoptim)

if (opts.meta):

    # check that summary statistics files are specified
    try:
        assert opts.statfiles is not None
    except AssertionError:
        print('Input error: You must specify --statinfo when using with --meta. See --help for details.')
        raise

    nstudies = len(opts.statfiles)

    try:
        assert nstudies > 0
    except AssertionError:
        print('Input error: Please specify the summary statistics file/s')
        raise

    statfiles = opts.statfiles
    featurefiles = opts.featurefiles


    if opts.params is not None:
        if opts.params[0] == 0.0:
            print("Input warning: pi cannot be equal to 0. Setting it to 0.000001")
            opts.params[0] = 0.000001
        mparams = hyperparameters.reverse(np.array(opts.params))
    else:
        mparams = None

    if opts.outdir is not None:
        outdir = os.path.realpath(opts.outdir)
    else:
        outdir = os.getcwd()
    if opts.file_prefix is not None:
        file_prefix = opts.file_prefix
    else:
        studyname = os.path.split(outdir)[1]
        file_prefix = '%s_meta' % studyname

    meta.optimize_hyperparameters(statfiles, featurefiles, outdir, file_prefix, opts.zmax, opts.muvar, mparams)
