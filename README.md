<img src="https://i.imgur.com/Pe9yoUO.png" alt="B-LORE logo" width="250"/>

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
1. Genotype files in [Oxford format](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html), for all loci of interest (e.g. Locus001.gen, Locus002.gen, etc.).
2. Sample file in [Oxford format](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html) (e.g. study1.sample)

For meta-analysis, it uses the following input:
1. Output files B-LORE summary statistics.
2. List of loci to be analyzed. This is a single file containing 2 columns with no header. The first column lists the name of the loci (e.g. Locus001, Locus002, etc.) and the second column is a binary number (1 or 0) indicating if it is a SNP locus (1) or a covariate locus (0). [Note: The summary statistics at each study outputs this file]
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
      
## Usage

### Quick start
- Clone the repository
- `cd example`
- `tar -zxvf input.tar.gz`
This will create an example input folder, with genotypes at 20 loci for 3 populations, a sample file for each population and ENCODE data for the 20 loci.
- `./commands.sh` to run B-LORE on the 3 populations to generate summary statistics, followed by a meta-analysis.

### Command line arguments
An executable file to run B-LORE is provided as `bin/blore`. This can used as follows:
```
blore [--help] [COMMAND] [OPTIONS]
```
There are 2 commands for B-LORE:
- `--summary` : for creating summary statistics of individual studies.
- `--meta` : for meta-analysis from summary statistics of multiple studies.

Each of these 2 commands takes different options, as described below.

#### blore --summary [OPTIONS]
Create summary statistics of individual studies. Valid options are:

Option | Description | Priority | Default value
:---   | :---        |:---      | :--
&#x2011;&#x2011;gen&nbsp;*filename(s)*  | Input genotype file(s), all loci should have separate genotype files and specified here (wildcards allowed) | Required    | --
&#x2011;&#x2011;sample&nbsp;*filename*  | Input sample file | Required | --
&#x2011;&#x2011;pheno&nbsp;*string*     | Name of the phenotype as it appears in the header of the sample file| Optional | `pheno`
&#x2011;&#x2011;regoptiom               | If specified, the variance of the regularizer will be optimized, otherwise it will be N(0, σ<sup>2</sup>) where σ is specified by `--reg` | Optional | --
&#x2011;&#x2011;reg&nbsp;*float*        | Value of the standard deviation (σ) of the regularizer | Optional | 0.01
&#x2011;&#x2011;pca&nbsp;*int*          | Number of principal components of the genotype to be included as covariates | Optional | 0
&#x2011;&#x2011;cov&nbsp;*string(s)*    | Name of covariate(s) as they appears in the header of the sample file, multiple covariates can be specified as space-separated strings | Optional | None
&#x2011;&#x2011;out&nbsp;*directory*    | Name of the output directory where summary statistics will be created | Optional | directory of the genotype files
&#x2011;&#x2011;prefix&nbsp;*string* | Prefix for the summary statistics files | Optional | `_summary`

#### blore --meta [OPTIONS]
Perform meta-analysis from summary statistics of multiple studies. Valid options are:

Option | Description | Priority | Default value
:---   | :---        |:---      | :--
&#x2011;&#x2011;input&nbsp;*filename*   | Input file containing list of loci to be analyzed together | Required | --
&#x2011;&#x2011;statdir&nbsp;*filename(s)*   | Input directory of B-LORE summary statistics | Required | --
&#x2011;&#x2011;feature&nbsp;*filename(s)*    | Input file(s) for genomic feature tracks | Optional      | --
&#x2011;&#x2011;params&nbsp;*floats* | Initial values of the hyperparameters, requires 4 space-separated floats corresponding to β<sub>π</sub> μ σ σ<sub>bg</sub>| Optional | 0.01 0.0 0.01 0.01
&#x2011;&#x2011;muvar | If specified, μ will be optimized, otherwise it will be fixed to the initial value (default 0) | Optional | --
&#x2011;&#x2011;zmax&nbsp;*int* | Maximum number of causal SNPs allowed | Optional | 2
&#x2011;&#x2011;out&nbsp;*directory*    | Name of the output directory where result files will be created | Optional | current directory
&#x2011;&#x2011;prefix&nbsp;*string* | Prefix for the meta-analysis output files | Optional | `_meta`

#### Example
- Clone the repository
- `cd example`
- `tar -zxvf input.tar.gz`
This will create an example input folder, with genotypes at 20 loci for 3 populations, a sample file for each population and ENCODE data for the 20 loci.

View `commands.sh` in your favorite editor to see the commands, and execute `./commands.sh` to run B-LORE on the 3 populations to generate summary statistics, followed by a meta-analysis.

## Citation
- Saikat Banerjee, Lingyao Zeng, Heribert Schunkert and Johannes Soeding (2017). [Bayesian multiple logistic regression for case-control GWAS.](https://doi.org/10.1101/198911 "B-LORE method details") bioRxiv.


## License
B-LORE is released under the GNU General Public License version 3. See LICENSE for more details. Copyright Johannes Soeding and Saikat Banerjee.

## Contact
[Saikat Banerjee](https://www.mpibpc.mpg.de/staff/42459)
