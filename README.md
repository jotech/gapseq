# gapseq
_Informed prediction and analysis of bacteria metabolic pathways and genome-scale networks_

_gapseq_ is designed to combine metabolic pathway analysis with metabolic network reconstruction and curation.
Based on genomic information and databases for pathways and reactions, _gapseq_ can be used for:
- prediction of metabolic pathways from various databases
- transporter inference
- metabolic model creation
- multi-step gap filling 


## Installation
### Ubuntu/Debian/Mint
```
sudo apt install ncbi-blast+ git libglpk-dev r-base-core exonerate bedtools barrnap
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

### Centos/Fedora/RHEL
```
sudo yum install ncbi-blast+ git glpk-devel BEDTools exonerate hmmer
git clone https://github.com/tseemann/barrnap.git
export PATH=$PATH:barrnap/bin/barrnap # needs to be permanent => .bashrc ?
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

### MacOS
using [homebrew](https://brew.sh)
```
brew install git glpk blast bedtools r brewsci/bio/barrnap
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

### Troubleshooting
- at least ncbi **blast** version 2.2.27 (4/2013) is needed. If your disribution only contains an older version, try to download a binary directly from [ncbi](https://shorturl.at/jkAH0)
- If you are getting the installation error ``'lib = "../R/library"' is not writable`` while installing the **R packages**, then try this command beforehand:
```
Rscript -e 'if( file.access(Sys.getenv("R_LIBS_USER"), mode=2) == -1 ) dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)'
```

## Quickstart
 This predicts network candidate reactions, builds a draft model and performs gap filling:
```
./gapseq doall toy/myb71.fna
```
Do the same but with a defined medium for gap filling:
```
./gapseq doall toy/ecoli.fna.gz dat/media/MM_glu.csv
```

## Pathway analysis
Search for chitin pathways:
```
./gapseq find -p chitin toy/myb71.fna.gz
```
Check for a certain enzyme availability, for example cytochrome c oxidases (oxidase test)
```
./gapseq find -e 1.9.3.1 toy/ecoli.fna.gz
```
Search for enzymes by name:
```
./gapseq find -r ligninase toy/myb71.fna.gz
```

## Creation and gap filling of metabolic models
1) Metabolic pathway analysis
```
./gapseq find -p all toy/myb71.fna
./gapseq find-transport toy/myb71.fna
```

2) Creation of draft model
```
 ./gapseq draft -r toy/myb71-all-Reactions.tbl -t toy/myb71-Transporter.tbl -p toy/myb71-all-Pathways.tbl -c toy/myb71.fna.gz
 ```

3) Gap filling
```
gapseq fill -m toy/myb71-draft.RDS -c toy/myb71-rxnWeights.RDS -n dat/media/TSBmed.csv
```
