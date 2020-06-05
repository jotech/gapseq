### Ubuntu/Debian/Mint
```
sudo apt install ncbi-blast+ git libglpk-dev r-base-core exonerate bedtools barrnap
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

### Centos/Fedora/RHEL
```
sudo yum install ncbi-blast+ git glpk-devel BEDTools exonerate hmmer
git clone https://github.com/tseemann/barrnap.git
export PATH=$PATH:barrnap/bin/barrnap # needs to be permanent => .bashrc ?
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

### MacOS
using [homebrew](https://brew.sh)
```
brew install git glpk blast bedtools r brewsci/bio/barrnap
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

### Troubleshooting
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
- we recommend using *cplex* als LP-solver. Please see [gapseq-Wiki](https://github.com/jotech/gapseq/wiki/Installation-of-cplexAPI) for installation instructions. 
