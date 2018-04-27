""" Optimize the mu and sigma of regularization prior used in quasi Laplace approximation
"""

import numpy as np
from scipy import optimize
from inference.logistic_regression import LogisticRegression
from inference.empirical_bayes import EmpiricalBayes

class OptimizeRegularizer:


    def __init__(self, genotype, phenotype, mu = 0.0, sigma = 0.01, tol = 0.01, covariates = None):
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
        sigmareg_old2 = self._sigmareg

        skip_step = False
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

            mureg_new = mu
            sigmareg_new = np.sqrt(sigma2)
            checkmu  = self.check_convergence(mureg_new, mureg_old)
            mureg_old = mureg_new

            if emp_bayes.success:

                # Check for convergence
                if (sigmareg_new > 0):
                    if skip_step:
                        checksig = self.check_convergence(sigmareg_new, sigmareg_old2)
                    else:
                        checksig = self.check_convergence(sigmareg_new, sigmareg_old)
                    sigmareg_old = sigmareg_new
                else:
                    # negative value of sigma 0. RESTART!
                    mureg_old = 0
                    sigmareg_old = np.random.normal(0.01, 0.001)

                skip_step = False # set it for the next round
    
                if checksig and checkmu:
                    iterate = False
                    mureg_optimized = mureg_new
                    sigmareg_optimized = sigmareg_new

            else:
                sigmareg_old2 = sigmareg_old
                sigmareg_old = sigmareg_new
                skip_step = True
                
            self._niter += 1

        self._mureg = mureg_optimized
        self._sigmareg = sigmareg_optimized


    def check_convergence(self, x, xold):
        check = False
        tol = self._tolerance
        diff = x - xold
        if diff == 0:
            check = True
        if not check and xold != 0.0:
            diff_percent = abs(diff) / xold
            if diff_percent < tol:
                check = True
            
        return check
