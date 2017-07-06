import numpy as np
from snptools.snpinfo import SnpInfo
import collections
from utils import matutils

class Editor:

    _standardize_once = False

    def __init__(self, genotype, snpinfo):
        self._genotype = list(genotype)
        self._snpinfo = list(snpinfo)
        self._nloci = len(genotype)
        self._low_maf = [[] for x in range(self._nloci)]
        self._duplicate = [[] for x in range(self._nloci)]


    @property
    def low_maf_snps(self):
        return tuple(self._low_maf)

    @property
    def duplicate_snps(self):
        return tuple(self._duplicate)

    @property
    def snpfreq(self):
        self._calculate_snpfreq()
        return tuple(self._snpfreq)

    @property
    def genotype(self):
        return tuple(self._genotype)

    @property
    def nsnps(self):
        return np.array([x.shape[1] for x in self._genotype])


    def remove_duplicates(self):
        for i in range(self._nloci):
            gt = self._genotype[i]
            nsnps = gt.shape[1]
            snps = self._snpinfo[i]

            unique_gt, unique_idx = matutils.unique_matrix(gt.T)
            newgt = unique_gt.T
            all_idx = set(np.arange(nsnps))
            removed_idx = all_idx - set(unique_idx)
            removed = list()
            for j in range(len(removed_idx)):
                removed.append(snps[list(removed_idx)[j]])

            self._genotype[i] = newgt
            self._duplicate[i]  = removed


    def remove_low_maf_snps(self):
        for i in range(self._nloci):
            gt = self._genotype[i]
            freq = np.sum(gt, axis=0) / (2 * gt.shape[0])
            nsnps = gt.shape[1]
            snps = self._snpinfo[i]

            removed = list()
            removed_idx = list()
            for j in range(nsnps):
                if (freq[j] > 0.99) or (freq[j] < 0.01):
                    removed_idx.append(j)
                    removed.append(snps[j])
            newidx = set(np.arange(nsnps)) - set(removed_idx)
            newgt = gt[:, list(newidx)]

            self._genotype[i] = newgt
            self._low_maf[i] = removed
    

    @staticmethod
    def norm_binom(gt, f):
        gt = (gt - (2 * f)) / np.sqrt(2 * f * (1 - f))
        return gt


    def _calculate_snpfreq(self):
        if (not self._standardize_once):
            self._snpfreq = list()
            for i in range(self._nloci):
                gt   = self._genotype[i]
                freq = np.sum(gt, axis=0) / (2 * gt.shape[0])
                self._snpfreq.append(freq)


    def standardize(self):
        ''' Once standardized, snpfreq calculation becomes difficult, and hence should be saved before standardizing.
            This function cannot be allowed to run more than once given an instance
        '''
        if self._standardize_once:
            return
        self._standardize_once = True

        self._snpfreq = list()
        for i in range(self._nloci):
            gt    = self._genotype[i]
            freq  = np.sum(gt, axis=0) / (2 * gt.shape[0])
            newgt = self.norm_binom(gt, freq)

            self._genotype[i] = newgt
            self._snpfreq.append(freq)
