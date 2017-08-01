import numpy as np

from snptools import dosageparser
from snptools import gtutils
from inference.logistic_regression import LogisticRegression
from inference.optimize_regularizer import OptimizeRegularizer
from iotools import io_covariates
from iotools.io_summary import WriteSummary
from utils import hyperparameters


def parse_dosage_data(samplefile, genotypefiles, phenotypename):
    parser = dosageparser.DosageParser()
    parser.set_samplefile(samplefile)
    parser.set_genotypefiles(genotypefiles)
    parser.set_phenotypename(phenotypename)
    data = parser.read()
    return data


def pcacomp(data, n):
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n, random_state = 0)
    pca.fit(data)
    comp = pca.transform(data)
    return comp


def summary(samplefile, genotypefiles, phenotypename, covariatenames, mureg, sigreg, tolerance, npca, outdir, file_prefix, regoptim):

    # Process the input data
    data = parse_dosage_data(samplefile, genotypefiles, phenotypename)
    nsample   = data.nsample   # number of samples
    nloci     = data.nloci     # number of loci
    snpinfo   = data.snpinfo   # list of namedtuples
    phenotype = data.phenotype # phenotype of samples, 1 for diseased, 0 for normal / numpy array
    genotype  = data.genotype  # list of genotype-matrix for all loci
    locnames  = data.locusnames

    # Remove the snps which are duplicates.
    gteditor = gtutils.Editor(genotype, snpinfo)
    gteditor.remove_duplicates()

    # Binomial standardization
    gteditor.standardize()

    # Get the final genotype
    gt = gteditor.genotype
    #nsnps = gteditor.nsnps #np.array([x.shape[1] for x in gt])

    # Covariates
    cov, covinfo = io_covariates.read(samplefile, covariatenames, nsample)
    if npca > 0:
        covinfo[0] += ["PC%i" % (i+1) for i in range(npca)]
        gt_tot = np.concatenate(gt, axis=1)
        pca = pcacomp(gt_tot, npca)
        pcastd = (pca - np.mean(pca, axis = 0)) / np.std(pca, axis=0)
        thiscov = pcastd
        if cov is not None:
            cov = np.append(cov, thiscov, axis=1) # cov will have a shape of nsample x ncov
        else:
            cov = thiscov

    # Optimize sigma_reg
    if regoptim:
        print ("Optimizing regularizer ...")
        cmax = 1
        sigreg_optim = OptimizeRegularizer(gt, phenotype, mureg, sigreg, tolerance, cov)
        sigreg_optim.update()
        mureg = sigreg_optim.mureg
        sigreg = sigreg_optim.sigmareg

    # Run the logistic regression
    print ("Calculating summary statistics")
    logreg = LogisticRegression(gt, phenotype, mureg, sigreg, covariates = cov)
    logreg.fit()
    v0 = logreg.v0
    vmin = logreg.vmin
    precll = logreg.precll
    iscov = logreg.iscov

    # Print out the results
    print ("Saving results")
    dupes = gteditor.duplicate_snps
    freq  = gteditor.snpfreq
    summary = WriteSummary(outdir, file_prefix)
    summary.set_statistics(snpinfo, covinfo, dupes, freq, v0, vmin, precll, mureg, sigreg, iscov, locnames)
    summary.write()
