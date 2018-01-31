import numpy as np
from scipy import optimize
from inference import zstatespy
from inference import quasi_laplace
from utils import hyperparameters

class EmpiricalBayes:
    ''' Empirical Bayes approach to learn the hyperparameters from summary statistics.
        Also called type-2 maximum likelihood or evidence approximation.
        We find the hyperparameters by maximising the marginal likelihood.

        Gets the marginal likelihood function and its derivatives from quasi_laplace.
        Uses a C implementation of creating z-states.
    '''


    _global_zstates = []
    _update_zstates = True
    _count_iter     = 0


    def __init__(self, vmin, precll, features, iscov, mureg, sigreg2, cmax, muvar, zcut = 0.98, params = None, regoptim = False, rerun = False):
        
        # Store all input variables
        self._vmin = vmin
        self._precll = precll
        self._features = features
        self._iscov = iscov
        self._mureg = mureg
        self._sigmareg2 = sigreg2
        self._cmax = cmax
        self._cutoff = zcut
        self._params = params if params is not None else np.array([-5.5, 0, -5, -5])
        self._rerun  = rerun
        self._regoptim = regoptim
        self._muvar = muvar

        # Create some more for internal use
        self._nloci = len(vmin)  # L + 1 (if covariates are present, else L)
        self._nsnps = [len(vmin[i]) for i in range(self._nloci)] # Number of variants in each locus
        self._nsnpfeat = [x.shape[0] for i, x in enumerate(features) if iscov[i] == False][0] # Get the number of features of each SNP locus, and use the first one
        self._ncovfeat = sum(iscov)

        # Different option for optimizing regularizer
        if regoptim:
            zstates = [[[i for i in range(x)]] for x in self._nsnps]
            self.fix_zstates(zstates)
            self._params = np.array([-5.5, 0, -7.5])


    @property
    def params(self):
        return self._params


    @property
    def zstates(self):
        return self._global_zstates


    def fit(self):

        # For the optimization of regularizer, all SNP features are 1, covariates are also considered as SNPs -- hence kc = 0
        k  = self._nsnpfeat if not self._regoptim else 1
        kc = max(1, self._ncovfeat) if not self._regoptim else 0

        # No need to create params if it is a rerun
        if self._rerun:
            params = self._params

        # Else, initialize the parameters
        else:
            params = np.zeros(k + 2 + kc)                     # k beta_pi for k features, 2 for beta_mu and beta_sigma, 1 for beta_cov
            params[0] = self._params[0]                       # baseline pi taken from input
            params[1 : k] = 0.0                               # beta_pi for all features are initialized to 0
            params[k] = self._params[1]                       # beta_mu taken from input
            params[k + 1] = self._params[2]                   # beta_sigm taken from input
            if kc > 0:
                params[-kc:] = self._params[3]                # beta_cov taken from input (only if kc > 0, note kc = 0 if regoptim)

        # Specify the bounds
        bounds = [[None, None] for i in range(k + 2 + kc)]
        if not self._muvar:
            for i in range(k, k + 1):
                bounds[i] = [params[i], params[i]]

        # Run the optimization
        args = self._nloci, self._vmin, self._precll, self._features, self._mureg, self._sigmareg2, self._iscov, self._regoptim, self._muvar
        self._callback_zstates(params)

        lml_min = optimize.minimize(self._log_marginal_likelihood,
                                    params,
                                    args = args,
                                    method='L-BFGS-B',
                                    jac=True,
                                    bounds=bounds,
                                    callback=self._callback_zstates,
                                    options={'maxiter': 200000,
                                             'maxfun': 2000000,
                                             'ftol': 1e-9,
                                             'gtol': 1e-9,
                                             'disp': True
                                            })

        self._params = lml_min.x
        #print (lml_min)


    def fix_zstates(self, zstates):
        self._global_zstates = zstates
        self._update_zstates = False


    def _log_marginal_likelihood(self, params, nloci, vmin, precll, feat, mureg, sigreg2, iscov, regoptim, muvar):

        zglobal = self._global_zstates
        nhyp = params.shape[0]  # Number of hyperparameters

        mll = 0
        der = np.zeros(nhyp)
        nfeat = feat[0].shape[0]

        der_list = list()
    
        for i in range(nloci):

            # Features for this locus
            # feat[i] => K x I for SNP loci, and Kc x J for covariate loci

            # Calculate pi, mu, sig2 for this locus from the hyperparameters and features
            pi, mu, sig2 = hyperparameters.transform(params, feat[i], iscov[i], regoptim)

            if regoptim:
                mll_locus, der_locus = self.reg_optim_mll_grad_locus(mu, sig2, vmin[i], nhyp, mureg, sigreg2, precll[i])
            else:
                mll_locus, der_locus = quasi_laplace.mll_grad(pi, mu, sig2, zglobal[i], vmin[i], nhyp, feat[i], mureg, sigreg2, precll[i], iscov[i])

    
            # Calculate the marginal log likelihood
            mll += mll_locus
            #print ("Locus MLL:", mll_locus)
    
            # Get the derivative of marginal log likelihood
            for j in range(nhyp):
                der[j] += der_locus[j]
            if not muvar:
                der[nfeat] = 0.0

            der_list.append(der_locus)

        #print (params)
        #print ("MLL: %g" % -mll)
        #print ("Gradient:", -der[2])

        # Priors for mu and beta_pi
        if not regoptim:
            blambda = 0.75
            logprior_pi = - blambda * np.sum(np.abs(params[1:nfeat]))
            derprior_pi = - blambda * np.sign(params[1:nfeat])
            #logprior_mu = - np.abs(params[nfeat])
            #derprior_mu = - np.sign(params[nfeat])
            musigma = 0.01
            mulambda = 1 / musigma / musigma
            logprior_mu = 0.5 * np.log(mulambda) - 0.5 * mulambda * params[nfeat] * params[nfeat] # Gaussian prior
            derprior_mu = - mulambda * params[nfeat]

            mll += logprior_pi
            mll += logprior_mu
            der[1:nfeat] += derprior_pi
            der[nfeat]  += derprior_mu

        return -mll, -der


    def _callback_zstates(self, params):
        if self._update_zstates and self._count_iter % 10 == 0:
            a = [hyperparameters.transform(params, self._features[i], self._iscov[i]) for i in range(self._nloci)]
            self._global_zstates = [zstatespy.create(self._nsnps[i],
                                                     self._cmax,
                                                     self._cutoff,
                                                     a[i][0], a[i][1], a[i][2], # pi, mu, sig2
                                                     self._vmin[i],
                                                     self._mureg,
                                                     self._sigmareg2,
                                                     self._precll[i],
                                                     self._iscov[i]) \
                                    for i in range(self._nloci)]
        self._count_iter += 1


    @staticmethod
    def reg_optim_mll_grad_locus(mu, sig2, lvmin, nhyp, mureg, sigreg2, lprecll):

        sigfact = 1 / sig2 - 1 / sigreg2
        nvar = len(lvmin)
        vzfact  = np.dot(lprecll, lvmin) + (mu / sig2) - (mureg / sigreg2)
        lamsz = np.array(lprecll)
        lamsz[np.diag_indices(nvar, ndim=2)] += sigfact
        lamsz_inv = np.linalg.inv(lamsz)
        musz = np.dot(lamsz_inv, vzfact)

        f1 = np.sum(mu * mu / sig2)
        f2 = np.dot(vzfact.T, musz)
        c1 = np.sum(np.log(1 / sig2))
        c2 = np.linalg.slogdet(lamsz)
        sign = c2[0]
        mll_locus = 0.5 * (c1 - c2[1] - f1 + f2)

        # Derivative wrapper
        der_locus = np.zeros(nhyp)
        der_locus[1] = np.sum((musz -  mu) / sig2)
        der_locus[2] = 0.5 * np.sum((np.square(musz - mu) - sig2 + np.diag(lamsz_inv)) / sig2)

        return mll_locus, der_locus
