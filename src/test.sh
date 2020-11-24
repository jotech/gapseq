#!/bin/bash

curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")

# temporary working directory
cd $(mktemp -d)

# gapseq version
echo -e `$dir/../gapseq -v`

# printing operation system
echo $OSTYPE
echo -e `uname -v` "\n"

echo "#######################"
echo "#Checking dependencies#"
echo "#######################"
i=0

check_cmd(){
    cmd=$1
    arg=$2
    print=$3
    if command -v $cmd &> /dev/null
    then
        version=$(eval $cmd "$arg")
    [[ ! "$print" = false ]] && echo $version
    else
        echo $cmd NOT FOUND
        i=$((i+1))
    fi
}

check_cmd awk "--version | head -n 1"
check_cmd sed "--version | head -n 1"
check_cmd grep "-V | head -n 1"
check_cmd perl "-v | head -n 2 | tail -n 1"
check_cmd tblastn "-version | head -n 1"
check_cmd exonerate "-v | head -n 1"
check_cmd bedtools "--version"
check_cmd barrnap "--help | head -n 2 | tail -n 1"
check_cmd R "--version | head -n 1"
check_cmd Rscript "--version | head -n 1" false
check_cmd git "--version"
echo -e "\nMissing dependencies: $i\n\n\n"

echo "#####################"
echo "#Checking R packages#"
echo "#####################"
Rscript -e 'needed.packages <- c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI", "BiocManager", "Biostrings", "jsonlite", "CHNOSZ"); avail.packages <- installed.packages(); i=0; for( pkg in needed.packages ){; idx <- match(pkg, avail.packages[,"Package"]); if( ! is.na(idx) ){; cat(pkg, avail.packages[idx,"Version"], "\n"); }else{; cat(pkg, "NOT FOUND", "\n"); i=i+1; }; }; cat("\nMissing R packages: ", i, "\n\n")'

echo "##############################"
echo "#Checking basic functionality#"
echo "##############################"
i=0

# Optimization
Rscript -e 'suppressMessages(library(sybil)); suppressMessages(library(glpkAPI)); data(Ec_core); sol <- sybil::optimizeProb(Ec_core); sol.ok <- sol@lp_stat == 5; cat("Optimization test:", ifelse( sol.ok, "OK", "FAILED"), "\n"); ' > opt.log
cat opt.log
if  grep -q "OK" opt.log; then
    i=$((i+1))
fi

# blast
fasta=$dir/../dat/seq/Bacteria/rev/1.2.4.1.fasta
tmp_file=$(mktemp)
gunzip -c $dir/../toy/myb71.fna.gz > $tmp_file
makeblastdb -in $tmp_file -dbtype nucl -out orgdb >/dev/null
tblastn -db orgdb -query $fasta -qcov_hsp_perc 75 -outfmt "6 $blast_format" > blast.out
if [[ -s blast.out ]]; then
    echo "Blast test: OK"
    i=$((i+1))
else
    echo "Blast test: FAILED"
fi

echo -e "\nPassed tests: $i/2\n"
