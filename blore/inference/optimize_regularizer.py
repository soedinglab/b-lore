""" Optimize the mu and sigma of regularization prior used in quasi Laplace approximation
"""

import numpy as np
from scipy import optimize
from inference.logistic_regression import LogisticRegression
from inference.empirical_bayes import EmpiricalBayes

class OptimizeRegularizer:


    def __init__(self, genotype, phenotype, mu = 0.0, sigma = 0.01, tol = 0.0001, covariates = None):
        self._genotype = genotype
        self._phenotype = phenotype
        self._mureg = mu
        self._sigmareg = sigma
        self._tolerance = tol
        self._cmax = 1
        self._covariates = covariates
        self._niter = 0


    @property
    def mureg(self):
        return self._mureg


    @property
    def sigmareg(self):
        return self._sigmareg


    @property
    def niter(self):
        return self._niter

    
    def update(self):
        mureg_old = self._mureg
        sigmareg_old = self._sigmareg
        tol = self._tolerance

        iterate = True
        while iterate:

            # Run the logistic regression with old regulariser
            logreg = LogisticRegression(self._genotype, self._phenotype, mureg_old, sigmareg_old, covariates=self._covariates)
            logreg.fit()
            v0 = logreg.v0
            vmin = logreg.vmin
            precll = logreg.precll
            iscov = logreg.iscov

            # Run the empirical Bayes with old regulariser and z = 1 [note ||z|| = nsnps]
            sigreg2_old = sigmareg_old * sigmareg_old
            features = [np.ones((1, len(x))) for x in vmin]
            muvar = False
            emp_bayes  = EmpiricalBayes(vmin, precll, features, iscov, mureg_old, sigreg2_old, self._cmax, muvar, regoptim = True)
            emp_bayes.fit()
            params = emp_bayes.params
            mu = params[1]
            sigma2 = np.exp(params[2])

            # Get new values for mu_reg and sigma_reg
            mureg_new = mu
            sigmareg_new = np.sqrt(sigma2)

            # Check for convergence
            mureg_diff = mureg_new - mureg_old
            sigmareg_diff = sigmareg_new - sigmareg_old
            #print ("New sigmareg = %g, old sigmareg = %g" % (sigmareg_new, sigmareg_old))
            #print ("**** Sigmareg difference is %g" % sigmareg_diff)
            #print ("New mureg = %g, old mureg = %g" % (mureg_new, mureg_old))
            #print ("**** Mureg difference is %g" % mureg_diff)
            if abs(mureg_diff) < tol and abs(sigmareg_diff) < tol :
                iterate = False
            mureg_old = mureg_new
            sigmareg_old = sigmareg_new

            self._niter += 1

        self._mureg = mureg_new
        self._sigmareg = sigmareg_new
