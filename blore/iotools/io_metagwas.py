import numpy as np
import os
from snptools.snpinfo import SnpInfo

class WriteResult:

    def __init__(self, outdir, file_prefix):
        self._outdir = outdir
        self._file_prefix = file_prefix
        self._statfile  = os.path.join(outdir, "%s.stat"       % file_prefix)
        self._logfile   = os.path.join(outdir, "%s.log"        % file_prefix)
        self._probdir   = os.path.join(outdir, "%s_res"        % file_prefix)
        self._hessfile  = os.path.join(outdir, "%s.prec.npy"   % file_prefix)
        self._covfile   = os.path.join(outdir, "%s.covariates" % file_prefix)
        self._zstatedir = os.path.join(outdir, "%s_zstates"    % file_prefix)

        if not os.path.isdir(self._outdir):
            os.makedirs(self._outdir)

        if not os.path.isdir(self._probdir):
            os.makedirs(self._probdir)

        if not os.path.isdir(self._zstatedir):
            os.makedirs(self._zstatedir)


    def set_result(self, resultlist, sigreg, mureg):
        self._resultlist = resultlist
        self._sigreg = sigreg
        self._mureg = mureg


    def write(self):
        nloci = 0
        ncov = 0
        nsnps = list()
        for res in self._resultlist:
            if not res.iscov:
                self.write_probfile(self._probdir, res)
                nloci += 1
                nsnps.append(len(res.finemap))
            else:
                self.write_covresult(self._covfile, res)
                ncov += 1
            self.write_zsfile(self._zstatedir, res)
        self.write_statfile(self._statfile, nloci, nsnps, ncov, self._mureg, self._sigreg)
        #np.save(self._hessfile, np.array(self._precll))


    def write_time(self, ttotal, tpreprocess, toptim, tout):
        with open(self._logfile, 'w') as mfile:
            mfile.write("Time taken\n")
            mfile.write("=============\n")
            mfile.write("\tTotal time: {:.2f} sec\n".format(ttotal))
            mfile.write("\tPreprocessing: {:2f} sec\n".format(tpreprocess))
            mfile.write("\tOptimization: {:.2f} sec\n".format(toptim))
            mfile.write("\tB-LORE time: {:2f} sec # Read data + Preprocessing + Optimization\n".format(toptim + tpreprocess + tread))
            mfile.write("\tWrite result: {:.2f} sec\n".format(tout))



    @staticmethod
    def write_statfile(statfile, nloci, nsnps, ncov, mureg, sigreg):
        nsnps_str = " ".join(["{:d}".format(i) for i in nsnps])
        with open(statfile, 'w') as mfile:
            mfile.write("{:d} loci\n".format(nloci))
            mfile.write("{:s}\n".format(nsnps_str))
            mfile.write("{:g}\t #sigreg: Variance of the Gaussian used as a regulariser for optimization\n".format(sigreg))
            mfile.write("{:g}\t #mureg: Mean of the Gaussian used as a regulariser for optimization\n".format(mureg))
            mfile.write("{:d} covariates\n".format(ncov))


    @staticmethod
    def write_probfile(outdir, res):
        outfile = os.path.join(outdir, "{:s}.res".format(res.locusname))
        filefmt = "{:15s} {:10s} {:s} {:s} {:g} {:g} {:g} {:g} {:g}\n"
        with open(outfile, 'w') as mfile:
            mfile.write("Causal Probability: {:g}\n".format(res.causal_prob))
            for snp in res.finemap:
                mfile.write(filefmt.format(snp.rsid, 
                                           snp.bp_location,
                                           snp.ref_allele,
                                           snp.alt_allele,
                                           snp.prob,
                                           snp.vmin,
                                           snp.pi,
                                           snp.mu,
                                           snp.sigma))


    @staticmethod
    def write_covresult(covfile, res):
        filefmt = "{:s} {:g} {:g} {:g} {:g} {:g}\n"
        with open(covfile, 'w') as mfile:
            mfile.write("Causal Probability: {:g}\n".format(res.causal_prob))
            for cov in res.finemap:
                mfile.write(filefmt.format(cov.rsid,
                                           cov.prob,
                                           cov.vmin,
                                           cov.pi,
                                           cov.mu,
                                           cov.sigma))


    @staticmethod
    def write_zsfile(outdir, res):
        outfile = os.path.join(outdir, "{:s}.zstates".format(res.locusname))
        with open(outfile, 'w') as mfile:
            for zres in res.zstates:
                zstr = ", ".join(["{:d}".format(x) for x in zres.index])
                mfile.write("[{:s}] {:g}\n".format(zstr, zres.zcomp))
