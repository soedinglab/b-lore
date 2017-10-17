# Parse data from dosage files

import numpy as np
from snptools.snpinfo import SnpInfo
import re
import os

class DosageParser:

    _read_phenotype_once = False
    _read_genotype_once = False
    _read_covariates_once = False
    _nloci = 0
    
    def set_samplefile(self, samplefile):
        self._samplefile = samplefile
        self._read_phenotype_once = False
        
    def set_genotypefiles(self, genotypefiles):
        self._list_of_genotypefiles = genotypefiles
        self._read_genotype_once = False
        self._nloci = len(self._list_of_genotypefiles)

    def set_phenotypename(self, phenoname):
        self._phenotypename = phenoname
        self._read_phenotype_once = False

    @property
    def nsample(self):
        self._read_phenotypes()
        return self._nsample
    
    @property
    def nloci(self):
        return self._nloci

    @property
    def snpinfo(self):
        self._run_once()
        return tuple(self._snpinfo)

    @property
    def nsnps(self):
        self._run_once()
        return tuple(self._nsnps)

    @property
    def genotype(self):
        self._read_genotypes()
        return tuple(self._genotype)

    @property
    def phenotype(self):
        self._read_phenotypes()
        return tuple(self._phenotype)

    @property
    def full_genotype(self):
        self._read_genotypes()
        return tuple(self._full_genotype)

    @property
    def age(self):
        self._read_covariates()
        return tuple(self._age)

    @property
    def sex(self):
        self._read_covariates()
        return tuple(self._sex)

    @property
    def locusnames(self):
        self._read_genotypes()
        return tuple(self._locusnames)

    def read(self):
        self._read_phenotypes()
        self._read_genotypes()
        return self


    def _read_phenotypes(self):
        if self._read_phenotype_once:
           return
        self._read_phenotype_once = True

        with open(self._samplefile, 'r') as samfile:
            header = samfile.readline().strip()
            header_strings = header.split()
            try: 
                ipheno = header_strings.index(self._phenotypename)
            except ValueError:
                print ('The specified phenotype %s is not in samplefile' % self._phenotypename)
                raise

        with open(self._samplefile, 'r') as samfile:
            phenotype = np.array([line.split()[ipheno] for line in samfile.readlines()[2:]], dtype=int)

        self._nsample = len(phenotype)
        self._phenotype = phenotype


    def _run_once(self):
        snpinfo = list()
        nsnps = list()
        locusnames = list()
        for i in range(self._nloci):
            snp_locus = list()
            filepath = self._list_of_genotypefiles[i]
            locusnames.append(os.path.basename(filepath))
            with open(filepath, 'r') as genfile:
                for line in genfile:
                    mline = line.split()
                    this_snp = SnpInfo(rsid = mline[1],
                                       bp_location = mline[2],
                                       alt_allele = mline[3],
                                       ref_allele = mline[4])
                    snp_locus.append(this_snp)

            nsnps_locus = len(snp_locus)
            snpinfo.append(snp_locus)
            nsnps.append(nsnps_locus)

        self._snpinfo = snpinfo
        self._nsnps = nsnps
        self._locusnames = locusnames
        nsample = (len(mline) - 5) / 3
        if float(nsample).is_integer():
            self._nsample = int(nsample)
        else:
            raise ValueError('Number of columns in genotype not divisible by 3')


    def _read_genotypes(self):
        if self._read_genotype_once:
            return
        self._read_genotype_once = True
        self._run_once() # otherwise, self._nsample and self._nsnps is not set
        genotype = list()
        dosage = list()
        for i in range(self._nloci):
            locus_dosage = self.read_single_locus_dosage(self._list_of_genotypefiles[i], self._nsample)
            locus_gt     = self.create_gt_matrix(locus_dosage, self._nsample, self._nsnps[i])
            dosage.append(locus_dosage)
            genotype.append(locus_gt)

        self._full_genotype = dosage
        self._genotype = genotype


    @staticmethod
    def read_single_locus_dosage(filename, nsample):
        print("Reading genotype from {:s}".format(os.path.basename(filename)))
        ncol = nsample * 3 + 5
        dosage = np.loadtxt(filename, usecols=range(5, ncol)).T
        return dosage


    @staticmethod
    def create_gt_matrix(dosage, nsample, nsnps):
        genotype = np.empty([nsample, nsnps], dtype=np.float64)
        for j in range(nsample):
            maj_ind = j * 3            # Allele AA, where A is the alternate allele
            het_ind = maj_ind + 1      # Allele AB, where B is the reference allele
            min_ind = het_ind + 1      # Allele BB
            genotype[j,:] = 2 * dosage[min_ind, :] + dosage[het_ind, :] # [AA, AB, BB] := [0, 1, 2]
        return genotype


    def _read_covariates(self):
        if self._read_covariates_once:
            return
        self._read_covariates_once = True
        abspath = os.path.abspath(self._samplefile)
        #infile = os.path.splitext(abspath)[0] + '.cov'
        infile = abspath
        age = [0 for x in range(self.nsample)]
        sex = [0 for x in range(self.nsample)]
        i = 0
        if os.path.exists(infile):
            with open(infile, 'r') as mfile:
                next(mfile)
                next(mfile)
                for line in mfile:
                    linesplit = line.split()
                    sex[i] = int(linesplit[5]) if not linesplit[5] == 'NA' else -1
                    age[i] = float(linesplit[6]) if not linesplit[6] == 'NA' else -1
                    i += 1
        self._age = np.array(age)
        self._sex = np.array(sex)

