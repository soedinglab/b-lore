#!/use/bin/env python

import numpy as np

def read(samplefile, covnames, nsample):

    covinfo = [[]]
    cov = None

    if covnames is not None:
        with open(samplefile, 'r') as samfile:
            header = samfile.readline().strip().split()
    
        for covname in covnames:
            try:
                icov = header.index(covname)
            except ValueError:
                print ('The specified covariate {:s} is not in samplefile'.format(covname))
            covinfo[0].append(covname) # currently we are using only one covariate locus.
            thiscov = [line.strip().split()[icov].strip() for line in open(samplefile, 'r').readlines()[2:]]
            A = [float(x) for x in thiscov if x != 'NA']
            mean = sum(A) / float(len(A))
            thiscov = [float(x) if x != 'NA' else mean for x in thiscov]
            thiscov = np.array(thiscov).reshape(nsample, 1)
            #covstd = (thiscov - np.mean(thiscov, axis = 0)) / np.std(thiscov, axis=0)
            covstd = thiscov
            if cov is not None:
                cov = np.append(cov, covstd, axis=1) # cov will have a shape of nsample x ncov
            else:
                cov = covstd

    return cov, covinfo
