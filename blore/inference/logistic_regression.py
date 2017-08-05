import numpy as np
from scipy import optimize

class LogisticRegression:


    _jac_calculated = False


    def __init__(self, genotype, phenotype, mureg, sigmareg, covariates=None):
        # list of genotype arrays each with shape (N,I)
        self._genotype = genotype
        # array of phenotypes shape(N,)
        self._phenotype = np.array(phenotype)
        # value of sigma regularizer
        self._sigmareg = sigmareg
        # value of mu regularizer
        self._mureg = mureg
        # array of covariates shape(N,J)
        self._covariates = covariates
        # Array for # of snps in each loci
        self._nsnps = np.array([x.shape[1] for x in genotype])


    @property
    def v0(self):
        return self._v0


    @property
    def vmin(self):
        return self._vmin


    @property
    def coeff(self):
        return self._coeff


    @property
    def pred(self):
        return self._predicted_risk


    @property
    def precll(self):
        self._jac()
        return self._precll


    @property
    def iscov(self):
        nloci = len(self._genotype)
        iscov = [False for i in range(nloci)]
        if self._covariates is not None:
            iscov += [True]
        return iscov


    def fit(self):

        # Prepare the input matrix for logistic regression
        # Concatenate genotype of all locus
        # Add covariates if any
        # Append a column of 1's for v0
        gtall = np.concatenate(self._genotype, axis=1)
        if self._covariates is not None:
            ncov_tot  = self._covariates.shape[1]
            gtcov = np.concatenate((gtall, self._covariates), axis=1)
        else:
            ncov_tot = 0
            gtcov = gtall
        ones  = np.ones(gtcov.shape[0]).reshape(gtcov.shape[0], 1)
        gtmod = np.concatenate((gtcov, ones), axis=1)

        phenotype = self._phenotype
        nsnps_tot = np.sum(self._nsnps)
        sigreg    = self._sigmareg
        #np.random.seed(0)
        vinit     = np.random.uniform(-0.001, 0.001, nsnps_tot + ncov_tot + 1)
        args      = gtmod, phenotype, sigreg

        gmode = optimize.minimize(self._log_likelihood,
                                  vinit,
                                  args=args,
                                  method='L-BFGS-B',
                                  jac=True,
                                  bounds=None,
                                  options={'maxiter': 20000000,
                                           'maxfun': 20000000,
                                           'ftol': 1e-9,
                                           'gtol': 1e-9
                                           #'disp': True
                                          })
        vmin_all = gmode.x
        vmin_snps = vmin_all[0:nsnps_tot]
        vmin_cov = vmin_all[nsnps_tot:nsnps_tot + ncov_tot]
        v0 = vmin_all[nsnps_tot + ncov_tot]
        vmin_split = np.split(vmin_snps, np.cumsum(self._nsnps))[:-1]

        # Separate vmin for each locus
        # Lists are mutable, so I prefer tuples
        # The covariates form a hypothetical locus
        self._coeff = vmin_all
        self._v0 = v0
        if len(vmin_cov) > 0:
            self._vmin = tuple(vmin_split + [vmin_cov])
        else:
            self._vmin = tuple(vmin_split)
        self._predicted_risk = self._transform(gtmod)

        # For every new fit, the jac needs to be calculated again
        self._jac_calculated = False


    @staticmethod
    def _log_likelihood(v, x, p, sigma_reg):

        dim = len(v)
        vx = np.einsum('i,ji->j', v, x)
        fact = 1 / sigma_reg / sigma_reg

        # Function
        llv = np.sum(np.einsum('i,i->i', p, vx) - np.log(1 + np.exp(vx)))
        reg = 0.5 * dim * np.log(fact) - 0.5 * fact * np.dot(v, v)
        loglik = llv + reg

        # Gradient
        pred = 1 / (1 + np.exp(-vx))
        der = np.einsum('i,ij->j',(p - pred), x) - (v * fact)

        return -loglik, -der


    def _transform(self, x):
        vx = np.einsum('i,ji->j', self.coeff, x)
        pred = 1 / (1 + np.exp(-vx))
        return pred


    def _jac(self):
        if self._jac_calculated:
            return
        self._jac_calculated = True

        gt = self._genotype
        if self._covariates is not None:
            gt = gt + (self._covariates, ) # hypothetical locus of covariates

        #self._precll = tuple([self._precision_matrix(v0, vmin[i], gt[i], sig) for i in range(nloci)])
        self._precll = tuple([self._precision_matrix(gt[i]) for i in range(len(gt))])


    def _precision_matrix(self, x):
        n = x.shape[1]
        hess = - np.einsum('i,i,ij,ik->jk', self.pred, (1 - self.pred), x, x)
        hess[np.diag_indices(n)] -= 1.0 / self._sigmareg / self._sigmareg
        return -hess
