import numpy as np

from snptools import dosageparser
from snptools import gtutils
from inference.logistic_regression import LogisticRegression
from inference.optimize_regularizer import OptimizeRegularizer
from inference import framingham
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
    pca = PCA(n_components=n)
    pca.fit(data)
    comp = pca.transform(data)
    return comp


def summary(samplefile, genotypefiles, phenotypename, mureg, sigreg, tolerance, npca, outdir, file_prefix, regoptim, add_fram):

    # Process the input data
    data = parse_dosage_data(samplefile, genotypefiles, phenotypename)
    nsample   = data.nsample   # number of samples
    nloci     = data.nloci     # number of loci
    snpinfo   = data.snpinfo   # list of namedtuples
    phenotype = data.phenotype # phenotype of samples, 1 for diseased, 0 for normal / numpy array
    genotype  = data.genotype  # list of genotype-matrix for all loci
    age       = data.age
    sex       = data.sex
    locnames  = data.locusnames

    print("Input files read.")

    # Remove the snps which are duplicates.
    gteditor = gtutils.Editor(genotype, snpinfo)
    gteditor.remove_duplicates()

    # Binomial standardization
    gteditor.standardize()

    # Get the final genotype
    gt = gteditor.genotype
    #nsnps = gteditor.nsnps #np.array([x.shape[1] for x in gt])

    # Covariates
    cov = None
    covinfo = [[]]
    covlocusnames = []
    if npca > 0:
        gt_tot = np.concatenate(gt, axis=1)
        covinfo = [["PC%i" % (i+1) for i in range(npca)]]
        pca = pcacomp(gt_tot, npca)
        pcastd = (pca - np.mean(pca, axis = 0)) / np.std(pca, axis=0)
        cov = pcastd
        print(cov)
        covlocusnames = ["{:d}pca_eigenvec".format(npca)]

    # Framingham risk
    if add_fram:
        fram = framingham.risk(nsample, age, sex)
    else:
        fram = np.zeros(nsample)

    # Optimize sigma_reg
    if regoptim:
        cmax = 1
        sigreg_optim = OptimizeRegularizer(gt, phenotype, fram, mureg, sigreg, tolerance, cov)
        sigreg_optim.update()
        mureg = sigreg_optim.mureg
        sigreg = sigreg_optim.sigmareg

    # Run the logistic regression
    logreg = LogisticRegression(gt, phenotype, mureg, sigreg, framingham = fram, covariates = cov)
    logreg.fit()
    v0 = logreg.v0
    vmin = logreg.vmin
    precll = logreg.precll
    iscov = logreg.iscov

    # Print out the results
    dupes = gteditor.duplicate_snps
    freq  = gteditor.snpfreq
    summary = WriteSummary(outdir, file_prefix)
    summary.set_statistics(snpinfo, covinfo, dupes, freq, v0, vmin, precll, mureg, sigreg, iscov, locnames, covlocusnames)
    summary.write()
