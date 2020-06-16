# Installation

## Ubuntu/Debian/Mint
```
sudo apt install ncbi-blast+ git libglpk-dev r-base-core exonerate bedtools barrnap
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

## Centos/Fedora/RHEL
```
sudo yum install ncbi-blast+ git glpk-devel BEDTools exonerate hmmer
git clone https://github.com/tseemann/barrnap.git
export PATH=$PATH:barrnap/bin/barrnap # needs to be permanent => .bashrc ?
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

## MacOS
using [homebrew](https://brew.sh)
```
brew install git glpk blast bedtools r brewsci/bio/barrnap
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

## Troubleshooting
- at least ncbi **blast** version 2.2.27 (4/2013) is needed. If your disribution only contains an older version, try to download a binary directly from [ncbi](https://shorturl.at/jkAH0)
- If you are getting the installation error ``'lib = "../R/library"' is not writable`` while installing the **R packages**, then try this command beforehand:
```
Rscript -e 'if( file.access(Sys.getenv("R_LIBS_USER"), mode=2) == -1 ) dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)'
```
- glpkAPI install error ``checking for library containing glp_create_prob... no`` although libglpk has been installed, try:
```
wget https://cran.r-project.org/src/contrib/glpkAPI_1.3.2.tar.gz
R CMD INSTALL --configure-args="--enable-gmp=no" glpkAPI_1.3.2.tar.gz
```
- we recommend using *cplex* als LP-solver as it is usually faster than *glpk*. The cplex solver is included in the *IBM ILOG CPLEX Optimization Studio*, which is free* for students and academics through the **IBM Academic Initiative** programm ([see here](https://developer.ibm.com/docloud/blog/2019/07/04/cplex-optimization-studio-for-students-and-academics/)). Please follow the installation instructions for *cplex* provided by IBM.
The R-package for the interface between R and the cplex solver can be optained from CRAN ([cplexAPI on CRAN](https://cran.r-project.org/web/packages/cplexAPI/index.html)).
Installation instructions (Linux & MacOS):
1. Download the package source archive from [CRAN](https://cran.r-project.org/web/packages/cplexAPI/index.html).
2. Open a terminal and **cd** to the directory of "cplexAPI_x.x.x.tar.gz"
3. Install die R-Package using the following commands:
```
CPLEXDIR=/dir/to/cplex129/cplex # adjust this with the path to your cplex installation
SYSTEM=`ls -d ${CPLEXDIR}/lib/*`
LIBFORMAT=static_pic
CPLEXLIBDIR=${SYSTEM}/${LIBFORMAT}

R CMD INSTALL --configure-args=" \
PKG_CFLAGS='-m64 -fPIC' \
PKG_CPPFLAGS='-I${CPLEXDIR}/include' \
PKG_LIBS='-L${CPLEXLIBDIR} -lcplex -lm -pthread -ldl'" \
cplexAPI_1.3.6.tar.gz
```

***
*source: [https://developer.ibm.com/docloud/blog/2019/07/04/cplex-optimization-studio-for-students-and-academics/](https://developer.ibm.com/docloud/blog/2019/07/04/cplex-optimization-studio-for-students-and-academics/)
