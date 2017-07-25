# B-LORE

Bayesian LOgistic REgression

A tool for meta-analysis in GWAS using Bayesian multiple logistic regression

## Description
B-LORE is a command line tool that creates summary statistics from multiple logistic regression on GWAS data,
and combines the summary statistics from multiple studies in a meta-analysis. 
It can also incorporate functional information about the SNPs from other external sources.
Several genetic regions, or loci are preselected for analysis with B-LORE.

### Key features
1. Association probability: B-LORE outputs probabilities of the input genetic loci being statistically associated with the phenotype.
2. Finemapping: B-LORE also outputs the probability of each SNP being statistically associated with the phenotype.
3. Leverage functional genomic data as a prior probability to improve prioritization.
4. Models data with logistic regression, and is suited for case/control studies.
5. Combines information over all SNPs in a locus with multiple regression.

## Installation
B-LORE is written in python and C++. To run B-LORE, you will need
- python version 3.4 or higher,
- the Python packages for scientific computing NumPy and SciPy.
- C++ compiler

To use B-LORE, you have to download the repository and compile the C++ shared libraries:
```
git clone https://github.com/soedinglab/b-lore.git
cd b-lore
make
```
The `Makefile` uses `g++` by default, which you can change depending on the compiler available on your system.


## Input files
For calculating summary statistics, it uses the following file formats as input:
1. Genotype files in [Oxford format](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html), for all loci of interest.
2. Sample file in [Oxford format](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html)

For meta-analysis, it uses the following input:
1. Output from B-LORE summary statistics. Note that it cannot use standard SNPTEST summary statistics.
2. (Optional) Functional genomics data, separately for each locus. 
Each feature file contains 2 parts: 
(a) a header line detailing the names of the columns in the file, and
(b) a line for each SNP detailing the information for that SNP.
The columns are tab-separated. 
The annotation tracks are present from column 4 onwards.
The first 3 columns are:
      - RSID: must have the same SNP identifier as in the genotype files
      - CHR: chromosome number
      - POS: base-pair position of the SNP.

## Quick start
- Clone the repository
- `cd example`
- `tar -zxvf input.tar.gz`
This will create an example input folder, with genotypes at 20 loci for 3 populations, a sample file for each population and ENCODE data for the 20 loci.
- `./commands.sh` to run B-LORE on the 3 populations to generate summary statistics, followed by a meta-analysis.

## License
B-LORE is released under the GNU General Public License version 3. See LICENSE for more details. Copyright Johannes Soeding and Saikat Banerjee.
