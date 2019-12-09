# gapseq
_Informed prediction and analysis of bacteria metabolic pathways and genome-scale networks_

_gapseq_ is designed to combine metabolic pathway analysis with metabolic network reconstruction and curation.
Based on genomic information and databases for pathways and reactions, _gapseq_ can be used for:
- prediction of metabolic pathways from various databases
- transporter inference
- metabolic model creation
- multi-step gap filling 


## Installation
```
sudo apt install ncbi-blast+ git libglpk-dev r-base-core exonerate bedtools barrnap
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq
cd gapseq
```

## Quickstart
 This performs the prediction of network candidate reactions, builds a draft model construction and performs gap filling.
```
./gapseq doall dat/myb71.fna
```

## Example workflow
1) Metabolic pathway analysis
```
./gapseq find -p all dat/myb71.fna
./gapseq find-transport.sh dat/myb71.fna
```

2) Creation of draft model
```
 ./gapseq draft -r toy/myb71-all-Reactions.tbl -t toy/myb71-Transporter.tbl -p toy/myb71-all-Pathways.tbl -c toy/myb71.fna.gz
 ```

3) Gap filling
```
gapseq fill -m toy/myb71-draft.RDS -c toy/myb71-rxnWeights.RDS -n dat/media/TSBmed.csv
```
