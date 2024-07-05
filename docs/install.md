# Installation

## Ubuntu/Debian/Mint
```
# Installation of main system dependencies
sudo apt install ncbi-blast+ git libglpk-dev r-base-core exonerate bedtools barrnap bc parallel curl libcurl4-openssl-dev libssl-dev libsbml5-dev bc

# installation of required R-packages
R -e 'install.packages(c("data.table", "stringr", "getopt", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI", "jsonlite", "httr"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'

# Download latest gapseq version from github
git clone https://github.com/jotech/gapseq && cd gapseq

# Download latest reference sequence database
bash ./src/update_sequences.sh
```

Test your installation with:
```sh
./gapseq test
```

## Centos/Fedora/RHEL
```
# Installation of main system dependencies
sudo yum install ncbi-blast+ git glpk-devel BEDTools exonerate hmmer bc parallel libcurl-devel curl openssl-devel libsbml-devel bc
git clone https://github.com/tseemann/barrnap.git
export PATH=$PATH:barrnap/bin/barrnap # needs to be permanent => .bashrc ?

# installation of required R-packages
R -e 'install.packages(c("data.table", "stringr", "getopt", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI", "CHNOSZ", "jsonlite", "httr"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
wget https://cran.r-project.org/src/contrib/Archive/sybil/sybil_2.2.0.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/sybilSBML/sybilSBML_3.1.2.tar.gz
R CMD INSTALL sybil_2.2.0.tar.gz
R CMD INSTALL sybilSBML_3.1.2.tar.gz
rm sybil_2.2.0.tar.gz
rm sybilSBML_3.1.2.tar.gz

# Download latest gapseq version from github
git clone https://github.com/jotech/gapseq && cd gapseq

# Download latest reference sequence database
bash ./src/update_sequences.sh
```

Test your installation with:
```sh
./gapseq test
```

## MacOS
Using [homebrew](https://brew.sh). Please note: Some Mac-Users reported difficulties to install gapseq on MacOS using the following commands. The issues are mainly due to some Mac-specific functioning of central programs such as sed, awk, and grep. If you are experiencing issues, we recommend to try to install gapseq in an own conda environment using the steps described [below](#conda).
```
# Installation of main system dependencies
brew install coreutils binutils git glpk blast bedtools r brewsci/bio/barrnap grep bc gzip parallel curl bc

# installation of required R-packages
R -e 'install.packages(c("data.table", "stringr", "getopt", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI", "CHNOSZ", "jsonlite", "httr"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
curl -O https://cran.r-project.org/src/contrib/Archive/sybil/sybil_2.2.0.tar.gz
R CMD INSTALL sybil_2.2.0.tar.gz
rm sybil_2.2.0.tar.gz

# Download latest gapseq version from github
git clone https://github.com/jotech/gapseq && cd gapseq

# Download latest reference sequence database
bash ./src/update_sequences.sh
```

Test your installation with:
```sh
./gapseq test
```

Some additional discussion and and trouble shooting can be found here: [1](https://apple.stackexchange.com/a/69332), [2](https://github.com/jotech/gapseq/issues/28), [3](https://github.com/jotech/gapseq/issues/143#issuecomment-1263349021).

## conda
There is no conda gapseq package available yet but all dependencies can be installed from conda channels without super user rights using the following steps:

**1. Install Mini-/Anaconda**

Follow the instructions provided by conda to Install [Anaconda/Miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

**2. Create conda environment for gapseq and adding package sources**

```sh
# Cloning the development version of gapseq
git clone https://github.com/jotech/gapseq
cd gapseq

# Create and activate a conda environment "gapseq-dev"
conda env create -n gapseq-dev --file gapseq_env.yml
conda activate gapseq-dev

# install one additional R-package
R -e 'install.packages("CHNOSZ", repos="http://cran.us.r-project.org")'

# Download & Install R-package 'sybilSBML'
wget https://cran.r-project.org/src/contrib/Archive/sybilSBML/sybilSBML_3.1.2.tar.gz
R CMD INSTALL --configure-args=" \
--with-sbml-include=$CONDA_PREFIX/include \
--with-sbml-lib=$CONDA_PREFIX/lib" sybilSBML_3.1.2.tar.gz
rm sybilSBML_3.1.2.tar.gz

# Download reference sequence data
bash ./src/update_sequences.sh
```

**3. Test the installation**

```sh
./gapseq test
```

## SBML support
The Systems Biology markup Language (SBML) can be used to exchange model files between gapseq and other programs.

The above installationn instructions for linux systems and using conda should already include the SBML support. If there were no errors during the installation, you are all set using gapseq with SBML format exports.

Occasionally, the installation can cause some issues that is why SBML is listed as optional dependency.
There should be a ``libsbml`` package (version 5.18.0 or later) available for most linux distributions:
```
sudo apt install libsbml5-dev # debian/ubuntu
sudo yum install libsbml-devel # fedora/centos
```
For MacOS, libsbml is not part of homebrew but an installation file can be downloaded from [here](https://sourceforge.net/projects/sbml/files/libsbml/5.18.0/stable/Mac%20OS%20X/). An instruction to install libSBML on an M1 Mac system has been suggested by a gapseq user [here](https://github.com/jotech/gapseq/issues/143#issuecomment-1263349021).

Please make sure, that `libsbml` is installed together with its [extensions "fbc" and "groups"](https://sbml.org/software/libsbml/).

Next, we need to install the R-package `sybilSBML`, version 3.1.2. Unfortunately, the package is currently not available anymore from the CRAN repository. However, the latest version can still be downloaded and installed from CRAN's archives:

```
wget https://cran.r-project.org/src/contrib/Archive/sybilSBML/sybilSBML_3.1.2.tar.gz
R CMD INSTALL sybilSBML_3.1.2.tar.gz
```
A common problem ist that the library path could not be found during installation of ``sybilSBML``. In this case, it may be necessary to specify the ``include`` and ``lib`` folder:
```
R CMD INSTALL --configure-args="--with-sbml-include=/path/to/libsbml-5.18.0/include/ --with-sbml-lib=/path/to/libsbml-5.18.0/lib/" sybilSBML_3.1.2.tar.gz
```
The sybilSBML archive is available at [gitlab](https://gitlab.cs.uni-duesseldorf.de/general/ccb/sybilSBML) together with more detailed installation information also at [gitlab](https://gitlab.cs.uni-duesseldorf.de/general/ccb/sybilSBML/-/blob/master/inst/INSTALL). Also, you may find more information on troubleshooting sybilSBML installation at this [gist-discussion](https://gist.github.com/dosorio/ea4baf66ee68821014d7dc6d92a48c55).

## cplex solver support

We recommend using *cplex* as LP-solver as it is usually faster than *glpk*. The cplex solver is included in the *IBM ILOG CPLEX Optimization Studio*, which might be at no charge to students and academics through the **IBM Academic Initiative** programm ([see here](https://community.ibm.com/community/user/ai-datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students)). Please follow the installation instructions for *cplex* provided by IBM.
The R-package for the interface between R and the cplex solver can be obtained from github ([SysBioChalmers/cplexAPI](https://github.com/SysBioChalmers/cplexAPI)). For *cplexAPI* installation please refer to instructions [here](https://github.com/SysBioChalmers/cplexAPI/blob/3070a6b60bb650919ee1b0db8b8223de99f88c3a/inst/INSTALL). Please note, that the R-package *cplexAPI* is not available from the CRAN repository anymore. We are currently working on a new CPLEX support solution that does not depend on *cplexAPI*.  


## Troubleshooting
- NCBI **blast** version 2.2.30 (10/2014) or **newer** is needed. If your disribution only contains an older version, try to download a binary directly from [ncbi](https://shorturl.at/jkAH0)
  * Older blast version could cause the ``Error: Unknown argument: "qcov_hsp_perc"``
- If you are getting the installation error ``'lib = "../R/library"' is not writable`` while installing the **R packages**, then try this command beforehand:
```
Rscript -e 'if( file.access(Sys.getenv("R_LIBS_USER"), mode=2) == -1 ) dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)'
```
- glpkAPI install error ``checking for library containing glp_create_prob... no`` although libglpk has been installed, try:
```
wget https://cran.r-project.org/src/contrib/glpkAPI_1.3.2.tar.gz
R CMD INSTALL --configure-args="--enable-gmp=no" glpkAPI_1.3.2.tar.gz
```
- *SBML-Export*: In order to export valid SBML files, it is required to have the R-Package [*sybilSBML*](https://cran.r-project.org/web/packages/sybilSBML/index.html) with version 3.1.2 or higher. In addtion libSBML-core-plus-packages (>= 5.16) with *'groups'* and *'fbc'* extensions are required. See also "SBML Support".

***
*source: [https://community.ibm.com/community/user/datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students](https://community.ibm.com/community/user/datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students)
