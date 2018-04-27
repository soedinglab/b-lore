#!/usr/bin/env python

import numpy as np
import collections
import os
from iotools.io_summary import ReadSummary

STUDYINFO_FIELDS = ['snpinfo', 'vmin', 'precll', 'sigreg']
class StudyInfo(collections.namedtuple('_StudyInfo', STUDYINFO_FIELDS)):
    __slots__ = ()


class SummaryStatistics:

    def __init__(self, studydirs, locusnamesfile):
        locusprefixes, iscov = self.read_locusnames(os.path.realpath(locusnamesfile))
        self._studydirs = studydirs
        self._locusprefixes = locusprefixes
        self._nloci = len(locusprefixes)
        self._nstudies = len(self._studydirs)
        self._iscov = iscov


    @property
    def snpinfo(self):
        return self._snpinfo
 

    @property
    def vmin(self):
        return self._vmin


    @property
    def precll(self):
        return self._precll


    @property
    def iscov(self):
        return self._iscov


    @property
    def sigreg(self):
        return self._sigreg


    @property
    def mureg(self):
        return self._mureg


    @property
    def locnames(self):
        return self._locusprefixes


    @property
    def covnames(self):
        return []
        #return [x.rsid for x in self._snpinfo[-1]]


    @staticmethod
    def read_locusnames(locusfile):
        locusprefixes = list()
        iscov = list()
        with open(locusfile, 'r') as mfile:
            for mlstr in mfile:
                mline = mlstr.split()
                locusprefix = mline[0]
                if len(mline) > 1:
                    flag = bool(int(mline[1]))
                else:
                    flag = True
                locusprefixes.append(locusprefix)
                iscov.append(not flag)
        return locusprefixes, iscov


    def combine(self):
        ''' Combine the study summaries of each locus one by one
        '''

        nstudies = self._nstudies
        snpinfo_list  = [[] for i in range(nstudies)]
        vmin_list     = [[] for i in range(nstudies)]
        precll_list   = [[] for i in range(nstudies)]
        sigreg_list   = [[] for i in range(nstudies)]

        sigfact = 0
        mufact  = 0

        # Read the summary statistics from all studies
        for i, studydir in enumerate(self._studydirs):
            print ('Reading summary statistics from {:s}'.format(studydir))
            summary = ReadSummary(studydir, self._locusprefixes)
            summary.read()

            snpinfo_list[i]  = summary.snpinfo
            vmin_list[i]     = summary.vmin
            precll_list[i]   = summary.precll
            sigreg_list[i]   = summary.sigreg
    
            sigfact += 1 / summary.sigreg / summary.sigreg
            mufact  += summary.mureg / summary.sigreg / summary.sigreg

        print ("Combining {:d} studies".format(nstudies))
        # Calculate some important statistics
        sigreg = np.sqrt(1 / sigfact)
        mureg = sigreg * sigreg * mufact
        nloci = self._nloci

        # Combine summary statistics of each locus
        snpinfo = list()
        precll = list()
        vmin = list()
        for l in range(nloci):
            print ('Locus %d / %d' % (l + 1, nloci))
            studies = list()
            for s in range(nstudies):
                study = StudyInfo(snpinfo = snpinfo_list[s][l], 
                                  vmin = vmin_list[s][l],
                                  precll = precll_list[s][l],
                                  sigreg = sigreg_list[s])
                studies.append(study)
            lsnpinfo, lprecll, lvmin = self.combine_single_locus(studies)
            snpinfo.append(lsnpinfo)
            precll.append(lprecll)
            vmin.append(lvmin)

        self._snpinfo = snpinfo
        self._precll = precll
        self._vmin = vmin
        self._sigreg = sigreg
        self._mureg = mureg


    def combine_single_locus(self, studies):
        ''' takes a list of objects as input:
            each object must have 4 elements: study.snpinfo, study.vmin, study.precll and study.sigreg
            where each element corresponds to a given locus.
   
            Returns combined vmin, precll and snpinfo for the given locus
        '''
        
        unique_rsidlist = list()
        unique_snplist  = list()
        vmin_dict   = collections.defaultdict(lambda:0)
        precll_dict = collections.defaultdict(lambda:0)

        for study in studies:
            unique_snps = [snp for snp in study.snpinfo if snp.rsid not in unique_rsidlist]
            unique_rsid = [snp.rsid for snp in unique_snps]
            unique_rsidlist += unique_rsid
            unique_snplist += unique_snps

            rsidlist = [snp.rsid for snp in study.snpinfo]
            lambda_vs = np.dot(study.precll, study.vmin)

            for i, rsid_i in enumerate(rsidlist):
                vmin_dict[rsid_i] += lambda_vs[i]
                for j, rsid_j in enumerate(rsidlist):
                    if j >= i:
                        precll_dict[rsid_i, rsid_j] += study.precll[i,j]

        for s, study in enumerate(studies):
            rsidlist = [snp.rsid for snp in study.snpinfo]
            for rsid_i in unique_rsidlist:
                if rsid_i not in rsidlist:
                    precll_dict[rsid_i, rsid_i] += 1.0 / (study.sigreg * study.sigreg)

        # precll needs to be converted from dictionary keys to index
        nsnps = len(unique_snplist)
        precll = np.zeros((nsnps, nsnps))
        for (rsid1, rsid2), val in precll_dict.items():
            i = unique_rsidlist.index(rsid1)
            j = unique_rsidlist.index(rsid2)
            precll[i,j] = val
            precll[j,i] = val

        # vmin needs to be converted from dictionary keys to index
        vmin = np.zeros(nsnps)
        for rsid, val in vmin_dict.items():
            i = unique_rsidlist.index(rsid)
            vmin[i] = val
        vmin = np.dot(np.linalg.inv(precll), vmin)

        return unique_snplist, precll, vmin
