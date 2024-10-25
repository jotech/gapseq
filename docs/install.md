# Installation

## Ubuntu/Debian/Mint
```
# Installation of main system dependencies
sudo apt install ncbi-blast+ git libglpk-dev r-base-core exonerate bedtools barrnap bc parallel curl libcurl4-openssl-dev libssl-dev libsbml5-dev bc

# installation of required R-packages
R -e 'install.packages(c("data.table", "stringr", "getopt", "R.utils", "stringi", "jsonlite", "httr", "pak"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
R -e 'pak::pkg_install("Waschina/cobrar")'

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
R -e 'install.packages(c("data.table", "stringr", "getopt", "R.utils", "stringi", "jsonlite", "httr", "pak"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
R -e 'pak::pkg_install("Waschina/cobrar")'

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
brew install coreutils binutils git glpk blast bedtools r brewsci/bio/barrnap grep bc gzip parallel curl bc brewsci/bio/libsbml

# installation of required R-packages
R -e 'install.packages(c("data.table", "stringr", "getopt", "R.utils", "stringi", "jsonlite", "httr"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
R -e 'pak::pkg_install("Waschina/cobrar")'

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

**Install Mini-/Anaconda**: Follow the instructions provided by conda to install [Anaconda/Miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Using conda, you can either install a specific release of gapseq or the latest development version of gapseq:

### Stable gapseq release using conda

Thanks to @cmkobel, a gapseq conda package is available for linux and osx platforms:
```sh
conda create -c conda-forge -c bioconda -n gapseq gapseq
```

### Development version using conda

The following commands create a conda environment for gapseq (named `gapseq-dev`) and installs gapseq along with all it's dependencies.

```sh
# Cloning the development version of gapseq
git clone https://github.com/jotech/gapseq
cd gapseq

# Create and activate a conda environment "gapseq-dev"
conda env create -n gapseq-dev --file gapseq_env.yml
conda activate gapseq-dev

# Download reference sequence data
bash ./src/update_sequences.sh
```

## SBML support
The Systems Biology markup Language (SBML) can be used to exchange model files between gapseq and other programs.

The above installation instructions for linux systems, MacOS, and using conda should already include the SBML support. If there were no errors during the installation, you should be all set using gapseq with SBML format exports.

Occasionally, the installation can cause some issues that is why SBML is listed as optional dependency.
There should be a ``libsbml`` package (version 5.18.0 or later) available for most linux distributions:
```
sudo apt install libsbml5-dev # debian/ubuntu
sudo yum install libsbml-devel # fedora/centos
```
For MacOS, libsbml is part of homebrew from [the `brewsci/bio` tap](https://github.com/brewsci/homebrew-bio/pkgs/container/bio%2Flibsbml).

If you want to manually install libsbml, please make sure, that `libsbml` is installed together with its [extensions "fbc" and "groups"](https://sbml.org/software/libsbml/). Those are extensions of the libsbml library, which cobrar requires.

## cplex solver support

We recommend using *cplex* as LP-solver as it is usually faster than *glpk*. The cplex solver is included in the *IBM ILOG CPLEX Optimization Studio*, which might be at no charge to students and academics through the **IBM Academic Initiative** program ([see here](https://community.ibm.com/community/user/ai-datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students)). Please follow the installation instructions for *cplex* provided by IBM.
The R-package for the interface between R and the cplex solver can be obtained from github ([Waschina/cobrarCPLEX](https://github.com/Waschina/cobrarCPLEX)). For *cobrarCPLEX* installation please refer to instructions [here](https://github.com/Waschina/cobrarCPLEX?tab=readme-ov-file#installation).

## Troubleshooting
- NCBI **blast** version 2.2.30 (10/2014) or **newer** is needed. If your distribution only contains an older version, try to download a binary directly from [ncbi](https://shorturl.at/jkAH0)
  * Older blast version could cause the ``Error: Unknown argument: "qcov_hsp_perc"``
- If you are getting the installation error ``'lib = "../R/library"' is not writable`` while installing the **R packages**, then try this command beforehand:
```
Rscript -e 'if( file.access(Sys.getenv("R_LIBS_USER"), mode=2) == -1 ) dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)'
```


***
*source: [https://community.ibm.com/community/user/datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students](https://community.ibm.com/community/user/datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students)
