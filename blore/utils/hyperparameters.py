import numpy as np

#def convert_to_readable(params):
#    pi     = 1 / (1 + np.exp(-params[0]))
#    mu     = params[1]
#    sigma2 = np.exp(params[2])
#    return pi, mu, sigma2

def reverse(target):
    params = np.zeros(target.shape[0])
    betapi = -np.log((1 / target[0]) - 1)
    betamu = target[1]
    betasig = 2 * np.log(target[2])
    if target.shape[0] > 3:
        betacov = 2 * np.log(target[3])
        params =  np.array([betapi, betamu, betasig, betacov])
    else:
        params = np.array([betapi, betamu, betasig])
    return params

    
def transform(params, features, is_covariate=False, all_causal=False):
    k = features.shape[0]
    if not is_covariate and not all_causal:
        betapi = params[0 : k]
        A = np.einsum('k, ki -> i', betapi, features)
        pi = 1 / (1 + np.exp(-A))
        if (pi == 0).any():
            kpi = np.where(pi == 0)
            pi[kpi] = 1e-200
            print ("Warning: Pi values equal to zero.")
        if (pi == 1).any():
            kpi = np.where(pi == 1)
            pi[kpi] = 0.999999
            print ("Warning: Pi values equal to one.")
        mu = np.full(features.shape[1], params[k])
        sig2 = np.full(features.shape[1], np.exp(params[k+1]))
    elif not is_covariate and all_causal:
        pi = np.ones(features.shape[1])
        mu = np.full(features.shape[1], params[k])
        sig2 = np.full(features.shape[1], np.exp(params[k+1]))
    elif is_covariate:
        beta_sigcov2 = params[-k:]
        nvar = features.shape[1]
        pi = np.ones(nvar)
        mu = np.zeros(nvar)
        sig2 = np.exp(np.dot(beta_sigcov2, features))
    return pi, mu, sig2
