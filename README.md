# gapseq
Informed prediction and analysis of bacteria metabolic pathways and genome-scale networks

## Installation
``
sudo apt install ncbi-blast+ git libglpk-dev r-base-core exonerate bedtools barrnap
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "Biostrings", "R.utils", "stringi"))'
git clone https://github.com/jotech/gapseq
cd gapseq
``

## Quickstart
``
./gapseq_quickstart.sh dat/myb71.fna
``

## Example workflow
1) Obtaining candidate reactions
``
./gapseq.sh -p -b 200 all dat/myb71.fna
./transporter.sh -b 200 dat/myb71.fna
``

2) Create draft model
``
Rscript src/generate_GSdraft.R myb71-all-Reactions.tbl -t myb71-Transporter.tbl dat/myb71.fna -u 200 -l 100 -a 2 -p myb71-Pathways.tbl
``

3) Gapfilling
``
Rscript gf.suite.R -m myb71.RDS -c myb71-rxnWeights.RDS -n dat/media/TSBmed.csv
``
