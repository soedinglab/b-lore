import numpy as np
import os
import sys
import glob
from snptools.snpinfo import SnpInfo

class WriteSummary:

    def __init__(self, outdir, file_prefix):
        self._statfile    = os.path.join(outdir, "{:s}.stat".format(file_prefix))
        self._summarydir  = os.path.join(outdir, "summary")
        self._covarfile   = os.path.join(outdir, "{:s}.cov".format(file_prefix))
        self._logfile     = os.path.join(outdir, "{:s}.log".format(file_prefix))
        self._namesfile   = os.path.join(outdir, "{:s}.locusnames".format(file_prefix))

        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        if not os.path.isdir(self._summarydir):
            os.makedirs(self._summarydir)


    def set_statistics(self, snpinfo, covinfo, dupes, freq, v0, vmin, precll, mureg, sigreg, iscov, locusnames, niter):
        self._snpinfo = snpinfo
        self._covinfo = covinfo
        self._dupes = dupes
        self._freq = freq
        self._v0 = v0
        self._vmin = vmin
        self._precll = precll
        self._sigreg = sigreg
        self._mureg = mureg
        self._is_covariate = iscov
        self._locusnames = locusnames
        self._niter = niter


    def write(self):
        nloci = len(self._snpinfo)
        ncovloci = np.sum(self._is_covariate)
        nsnps = [len(self._snpinfo[i]) for i in range(nloci)]
        nsnps_str = " ".join(["{:d}".format(i) for i in nsnps])
        with open(self._statfile, 'w') as mfile:
            mfile.write("{:d} loci\n".format(nloci))
            mfile.write("{:s}\n".format(nsnps_str))
            mfile.write("{:g} #sigreg\n".format(self._sigreg))
            mfile.write("{:g} #mureg\n".format(self._mureg))
            mfile.write("{:g} #base_v0\n".format(self._v0))
            mfile.write("{:d} #Number of covariate loci\n".format(ncovloci))
            mfile.write("{:d} #Number of iterations for optimizing regularizer\n".format(self._niter))

        for i, locus in enumerate(self._locusnames):
            summaryfile = os.path.join(self._summarydir, '{:s}.vmin'.format(locus))
            with open(summaryfile, 'w') as mfile:
                mfile.write("#RSID\tPOS\tALT\tREF\tDUPE\tFREQ_REF\tVMIN\n")
                snps = self._snpinfo[i]
                removed = [item.rsid for item in self._dupes[i]]
                k = 0
                for snp in snps:
                    if snp.rsid in removed:
                        indx = removed.index(snp.rsid)
                        mfile.write("{:15s}\t{:10s}\t{:s}\t{:s}\t{:1d}\tNA\tNA\n".format(snp.rsid, snp.bp_location, snp.alt_allele, snp.ref_allele, 1))
                    else:
                        mfile.write("{:15s}\t{:10s}\t{:s}\t{:s}\t{:1d}\t{:g}\t{:g}\n".format(snp.rsid, 
                                                                                             snp.bp_location,
                                                                                             snp.alt_allele,
                                                                                             snp.ref_allele,
                                                                                             0,
                                                                                             self._freq[i][k],
                                                                                             self._vmin[i][k]))
                        k += 1
            hessfile = os.path.join(self._summarydir, '{:s}.jac'.format(locus))
            np.savetxt(hessfile, self._precll[i])

        for i, covs in enumerate(self._covinfo):
            if len(covs) > 0:
                summaryfile = os.path.join(self._summarydir, 'Covariates_{:d}.vmin'.format(i+1))
                with open(summaryfile, 'w') as mfile:
                    mfile.write("#RSID\tPOS\tALT\tREF\tDUPE\tFREQ_REF\tVMIN\n")
                    for j, cov in enumerate(covs):
                        mfile.write("{:s}\t0\tNA\tNA\t0\t0\t{:g}\n".format(cov, self._vmin[nloci+i][j]))
                hessfile = os.path.join(self._summarydir, 'Covariates_{:d}.jac'.format(i+1))
                np.savetxt(hessfile, self._precll[nloci+i])

        with open(self._namesfile, 'w') as mfile:
            for locus in self._locusnames:
                mfile.write("{:s}\t1\n".format(locus))
            for i, covs in enumerate(self._covinfo):
                if len(covs) > 0:
                    mfile.write("Covariates_{:d}\t0".format(i+1))


    def write_time(self, ttotal, tread, tpreprocess, toptim, tlogreg, tout):
        with open(self._logfile, 'w') as mfile:
            mfile.write("Time taken\n")
            mfile.write("=============\n")
            mfile.write("\tTotal time: {:.2f} sec\n".format(ttotal))
            mfile.write("\tRead data: {:.2f} sec\n".format(tread))
            mfile.write("\tPreprocessing: {:2f} sec\n".format(tpreprocess))
            mfile.write("\tOptimization: {:.2f} sec\n".format(toptim))
            mfile.write("\tLogistic regression: {:2f} sec\n".format(tlogreg))
            mfile.write("\tB-LORE time: {:2f} sec # Optimization + Logistic regression\n".format(toptim + tlogreg))
            mfile.write("\tWrite result: {:.2f} sec\n".format(tout))


class ReadSummary:

    _read_vmin_once = False
    _read_stats_once = False

    def __init__(self, outdir, locusprefixes):

        thisdir = os.path.realpath(outdir)
        statfile = glob.glob(thisdir + '/*.stat')
        if len(statfile) == 0:
            print("No stat file found in {:s}".format(thisdir))
            sys.exit()
        elif len(statfile) > 1:
            print("Multiple stat files found in {:s}".format(thisdir))
            sys.exit()
        else:
            self._statfile = os.path.realpath(statfile[0])
        self._summarydir = os.path.join(outdir, "summary")
        self._locusprefixes = locusprefixes


    @property
    def snpinfo(self):
        self._read_vmin()
        return self._snpinfo


    @property
    def vmin(self):
        self._read_vmin()
        return self._vmin


    @property
    def precll(self):
        self._read_vmin()
        return self._precll


    @property
    def sigreg(self):
        self._read_stats()
        return self._sigreg


    @property
    def mureg(self):
        self._read_stats()
        return self._mureg


    @property
    def snpfreq(self):
        self._read_vmin()
        return self._freq


    @property
    def nsnps(self):
        self._read_stats()
        return self._nsnps


    @property
    def v0(self):
        self._read_stats()
        return self._v0


    def read(self):
        self._read_stats()
        self._read_vmin()


    def _read_stats(self):

        if self._read_stats_once:
            return
        self._read_stats_once = True

        try:
            with open(self._statfile, 'r') as mfile:
                next(mfile)
                nsnps_str    = mfile.readline().split()
                self._nsnps  = [int(i) for i in nsnps_str]
                self._sigreg = float(mfile.readline().split()[0])
                self._mureg  = float(mfile.readline().split()[0])
                self._v0     = float(mfile.readline().split()[0])
        except OSError as err:
            print ("OS error: {0}".format(err))
        except ValueError:
            print ("Error in input values in %s" % self._statfile)
        except:
            print ("Unexpected error while reading summary statistics from %s" % self._statfile)
            raise

    def _read_vmin(self):

        if self._read_vmin_once:
            return
        self._read_vmin_once = True

        self._read_stats()

        nloci = len(self._locusprefixes)
        self._vmin = [[] for i in range(nloci)]
        self._snpinfo = [[] for i in range(nloci)]
        self._freq = [[] for i in range(nloci)]
        jaclist = [None for i in range(nloci)]

        for l, locus in enumerate(self._locusprefixes):
            summaryfile = os.path.join(self._summarydir, '{:s}.vmin'.format(locus))
            hessfile    = os.path.join(self._summarydir, '{:s}.jac'.format(locus))
            try:
                with open(summaryfile, 'r') as mfile:
                    next(mfile)
                    for mlstr in mfile:
                        mline = mlstr.split()
                        if int(mline[4]) == 1: continue
                        this_snp = SnpInfo(rsid = mline[0],
                                           bp_location = mline[1],
                                           ref_allele = mline[3],
                                           alt_allele = mline[2])
                        self._snpinfo[l].append(this_snp)
                        self._freq[l].append(float(mline[5]))
                        self._vmin[l].append(float(mline[6]))
            except Exception as e:
                print ("Unexpected error while reading summary statistics from %s" % self._summaryfile)
                print (e)
            jac = np.loadtxt(hessfile, ndmin = 1)
            if jac.shape[0] == 1:
               jac = jac.reshape(1,1)
            jaclist[l] = jac

        self._vmin = tuple(self._vmin)
        self._snpinfo = tuple(self._snpinfo)
        self._freq = tuple(self._freq)
        self._precll = tuple(jaclist)
