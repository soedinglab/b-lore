""" Optimize the mu and sigma of regularization prior used in quasi Laplace approximation
"""

import numpy as np
from scipy import optimize
from inference.logistic_regression import LogisticRegression
from inference.empirical_bayes import EmpiricalBayes

class OptimizeRegularizer:


    def __init__(self, genotype, phenotype, framingham, mu = 0.0, sigma = 0.01, tol = 0.0001, covariates = None):
        self._genotype = genotype
        self._phenotype = phenotype
        self._framingham = framingham
        self._mureg = mu
        self._sigmareg = sigma
        self._tolerance = tol
        self._cmax = 1
        self._covariates = covariates


    @property
    def mureg(self):
        return self._mureg


    @property
    def sigmareg(self):
        return self._sigmareg

    
    def update(self):
        mureg_old = self._mureg
        sigmareg_old = self._sigmareg
        tol = self._tolerance

        iterate = True
        while iterate:

            # Run the logistic regression with old regulariser
            logreg = LogisticRegression(self._genotype, self._phenotype, mureg_old, sigmareg_old, framingham = self._framingham, covariates=self._covariates)
            logreg.fit()
            v0 = logreg.v0
            vmin = logreg.vmin
            precll = logreg.precll
            iscov = logreg.iscov

            # Run the empirical Bayes with old regulariser and z = 1 [note ||z|| = nsnps]
            sigreg2_old = sigmareg_old * sigmareg_old
            k = 1
            kc = int(np.sum(iscov))
            features = [np.ones((1, len(x))) for x in vmin]
            emp_bayes = EmpiricalBayes(vmin, precll, features, iscov, mureg_old, sigreg2_old, self._cmax, False, regoptim = True)
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
            print ("New sigmareg = %g, old sigmareg = %g" % (sigmareg_new, sigmareg_old))
            print ("**** Sigmareg difference is %g" % sigmareg_diff)
            if abs(mureg_diff) < tol and abs(sigmareg_diff) < tol :
                iterate = False
            mureg_old = mureg_new
            sigmareg_old = sigmareg_new

        self._mureg = mureg_new
        self._sigmareg = sigmareg_new
