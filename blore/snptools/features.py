#!/usr/bin/env python

import numpy as np
import collections

class FunctionalAnnotations:

    _read_features_once = False

    def __init__(self, snpinfo, iscov, filepaths):
        self._snpinfo = snpinfo
        self._iscov = iscov
        self._filepaths = filepaths
        self._nsnps = [len(x) for x in snpinfo]
        self._nloci = len(snpinfo)


    @property
    def features(self):
        self.read()
        return self._features


    def read(self):

        if self._read_features_once:
            return
        self._read_features_once = True

        # read the features if filepaths are provided, else return a vector of ones.
        if self._filepaths is None:
            featurelist = [np.ones((1, len(x))) for x in self._snpinfo]
        else:
            featurelist = list()
            for l, locusinfo in enumerate(self._snpinfo):
                # if this is a covariate locus, set all features to baseline
                if self._iscov[l]:
                    features = np.ones((1, self._nsnps[l]))
                else:
                    with open(self._filepaths[l], 'r') as featfile:
                        featnames = featfile.readline().split()[3:]
                        nfeat = len(featnames) + 1
                        features = np.zeros((nfeat, len(locusinfo)))  # initialize the array
                        features[0, :] = 1 # baseline feature of 1

                        # first store all features as vectors in dictionaries for a couple of checks (see below)
                        featdict = collections.defaultdict(lambda:0)
                        rsidlist = list()
                        bplist = list()
                        for line in featfile:
                            linesplit = line.split()
                            rsid = linesplit[0]
                            bp = linesplit[2]
                            rsidlist.append(rsid)
                            bplist.append(bp)
                            featdict[rsid] = np.array(linesplit[3:])

                        # once reading is complete, populate the features array
                        for i, snp in enumerate(locusinfo):
                            # check number of times this snp is present in feature file
                            ncount = rsidlist.count(snp.rsid)
                            if ncount == 0:
                                print ('Warning: SNP %s is absent in the feature file' % snp.rsid)
                            if ncount > 1:
                                print ('Warning: SNP %s is present multiple times in the feature file' % snp.rsid)
                            if ncount == 1:
                                features[1:, i] = featdict[snp.rsid]

                        # scale the features
                        #for i in range(1, features.shape[0]):
                        #    x = features[i, :]
                        #    x = (x - np.mean(x)) / np.std(x)
                        #    features[i, :] = x

                # Reading this locus is complete, append it to featurelist
                #features = features[:2, :]
                featurelist.append(features)

        self._features = featurelist
