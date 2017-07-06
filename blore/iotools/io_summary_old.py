import numpy as np
import os
from snptools.snpinfo import SnpInfo

class WriteSummary:

    def __init__(self, outdir, file_prefix):
        self._statfile    = os.path.join(outdir, "%s.stat"    % file_prefix)
        self._summaryfile = os.path.join(outdir, "%s.vmin"    % file_prefix)
        self._covarfile   = os.path.join(outdir, "%s.cov"     % file_prefix)
        self._hessfile    = os.path.join(outdir, "%s.jac.npy" % file_prefix)

        if not os.path.isdir(outdir):
            os.makedirs(outdir)


    def set_statistics(self, snpinfo, covinfo, dupes, freq, v0, vmin, precll, mureg, sigreg, iscov):
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


    def write(self):
        nloci = len(self._snpinfo)
        ncovloci = np.sum(self._is_covariate)
        nsnps = [len(self._snpinfo[i]) for i in range(nloci)]
        nsnps_str = ["%i" % i for i in nsnps]
        with open(self._statfile, 'w') as mfile:
            mfile.write("%i loci\n" % nloci)
            mfile.write(" ".join(nsnps_str) + '\n')
            mfile.write("%g #sigreg\n" % self._sigreg)
            mfile.write("%g #mureg\n" % self._mureg)
            mfile.write("%g #base_v0\n" % self._v0)
            mfile.write("%g #Number of hypothetical covariate loci\n" % ncovloci)

        with open(self._summaryfile, 'w') as mfile:
            mfile.write('Locus rsid bp_location alt_all ref_all is_dupe freq_ref vmin\n')
            for i in range(nloci):
                #mfile.write("%s" % self._locus_names[i])
                snps = self._snpinfo[i]
                removed = [item.rsid for item in self._dupes[i]]
                k = 0
                for j in range(len(snps)):
                    if snps[j].rsid in removed:
                        indx = removed.index(snps[j].rsid)
                        mfile.write("%3i\t%15s\t%10s\t%s\t%s\t%1i\tNA\tNA\n" % (i+1, snps[j].rsid, snps[j].bp_location, snps[j].alt_allele, snps[j].ref_allele, 1))
                    else:
                        mfile.write("%3i\t%15s\t%10s\t%s\t%s\t%1i\t%g\t%g\n" % (i+1, snps[j].rsid, snps[j].bp_location, snps[j].alt_allele, snps[j].ref_allele, 0, self._freq[i][k], self._vmin[i][k]))
                        k += 1

        if ncovloci > 0:
            with open(self._covarfile, 'w') as mfile:
                mfile.write('Locus Covariate_name vmin\n')
                for i in range(ncovloci):
                    cov = self._covinfo[i]
                    for j in range(len(cov)):
                        mfile.write("%3i\t%s\t%g\n" %(nloci+i+1, cov[j], self._vmin[nloci+i][j]))

        np.save(self._hessfile, np.array(self._precll))


class ReadSummary:

    _read_vmin_once = False
    _read_stats_once = False
    _read_hessian_once = False

    def __init__(self, infile):
        self._statfile    = os.path.realpath("%s.stat"    % infile)
        self._summaryfile = os.path.realpath("%s.vmin"    % infile)
        self._covarfile   = os.path.realpath("%s.cov"     % infile)
        self._hessfile    = os.path.realpath("%s.jac.npy" % infile)


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
        self._read_hessian()
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
    def iscov(self):
        self._read_vmin()
        return self._iscov


    @property
    def snpfreq(self):
        self._read_vmin()
        return self._freq


    @property
    def nloci(self):
        self._read_stats()
        return self._nloci + self._ncovloci


    @property
    def nsnps(self):
        self._read_stats()
        return self._nsnps


    @property
    def v0(self):
        self._read_stats()
        return self._v0


    @property
    def locnames(self):
        self._read_vmin()
        return self._locnames


    @property
    def covnames(self):
        self._read_vmin()
        covnames = list()
        if self._ncovloci > 0:
            covlist = [x for i, x in enumerate(self._snpinfo) if self._iscov[i]]
            for cov in covlist:
                covnames.append([x.rsid for x in cov])
        return covnames


    def read(self):
        self._read_stats()
        self._read_vmin()
        self._read_hessian()


    def _read_stats(self):

        if self._read_stats_once:
            return
        self._read_stats_once = True

        try:
            with open(self._statfile, 'r') as mfile:
                self._nloci  = int(mfile.readline().split()[0])
                nsnps_str    = mfile.readline().split()
                self._nsnps  = [int(i) for i in nsnps_str]
                self._sigreg = float(mfile.readline().split()[0])
                self._mureg  = float(mfile.readline().split()[0])
                self._v0     = float(mfile.readline().split()[0])
                self._ncovloci = int(mfile.readline().split()[0])
        except OSError as err:
            print ("OS error: {0}".format(err))
        except ValueError:
            print ("Error in input values in %s" % self._statfile)
        except:
            print ("Unexpected error while reading summary statistics from %s" % self._statfile)
            raise
        # Check the summary values
        if not self._nloci == len(self._nsnps):
            raise ValueError("Number of loci reported in first line is different from the number of list elements in second line of %s" % self._statfile)

    def _read_vmin(self):

        if self._read_vmin_once:
            return
        self._read_vmin_once = True

        self._read_stats()

        self._vmin = [[] for i in range(self._nloci + self._ncovloci)]
        self._snpinfo = [[] for i in range(self._nloci + self._ncovloci)]
        self._freq = [[] for i in range(self._nloci + self._ncovloci)]
        self._iscov = [False for i in range(self._nloci)] + [True for i in range(self._ncovloci)]
        self._locnames = ["Locus.{:03d}".format(i + 1) for i in range(self._nloci + self._ncovloci)]

        try:
            with open(self._summaryfile, 'r') as mfile:
                next(mfile)
                for line in mfile:
                    mline = line.split()
                    if int(mline[5]) == 0:
                        locus_indx = int(mline[0]) - 1
                        this_snp = SnpInfo(rsid = mline[1],
                                       bp_location = mline[2],
                                       alt_allele = mline[3],
                                       ref_allele = mline[4])
                        self._snpinfo[locus_indx].append(this_snp)
                        self._freq[locus_indx].append(float(mline[6]))
                        self._vmin[locus_indx].append(float(mline[7]))
        except:
            print ("Unexpected error while reading summary statistics from %s" % self._summaryfile)

        #self._covnames = [[] for i in range(self._ncovloci)]
        if self._ncovloci > 0:
            try:
                with open(self._covarfile, 'r') as mfile:
                    next(mfile)
                    for line in mfile:
                        mline = line.split()
                        locus_indx = int(mline[0]) - 1
                        this_cov = SnpInfo(rsid = mline[1],
                                           bp_location = "0",
                                           alt_allele = "NA",
                                           ref_allele = "NA")
                        #cov_indx = locus_indx - self._nloci
                        #self._covnames[cov_indx].append(mline[1])
                        self._snpinfo[locus_indx].append(this_cov)
                        self._freq[locus_indx].append(0.0)
                        self._vmin[locus_indx].append(float(mline[2]))
            except:
                print("Unexpected error while reading covariance statistics from %s" % self._self._covarfile)

        self._vmin = tuple(self._vmin)
        self._snpinfo = tuple(self._snpinfo)
        self._freq = tuple(self._freq)


    def _read_hessian(self):
        
        if self._read_hessian_once:
            return
        self._read_hessian_once = True

        self._precll = np.load(self._hessfile)
