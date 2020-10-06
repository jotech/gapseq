#!/bin/bash

curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")

# gapseq version
echo -e `$dir/../gapseq -v` "\n"

echo "#######################"
echo "#Checking dependencies#"
echo "#######################"
i=0

check_cmd(){
    cmd=$1
    arg=$2
    if [[ "$3" == "F" ]]; then # add no space
        space=""
    else
        space="  "
    fi
    #space=$3
    #[[ -z "$space" ]] && { space="  "; echo empty; }
    if command -v $cmd &> /dev/null
    then
        version=$(eval $cmd "$arg")
    echo "$space $version"
    else
        echo "$space $cmd NOT FOUND"
        i=$((i+1))
    fi
}

check_cmd awk "-V | head -n 1"
check_cmd sed "--version | head -n 1"
check_cmd grep "-V | head -n 1"
check_cmd perl "-v | head -n 2 | tail -n 1"
check_cmd tblastn "-version | head -n 1"
check_cmd exonerate "-v | head -n 1"
check_cmd bedtools "--version"
check_cmd barrnap "--help | head -n 2 | tail -n 1" F
check_cmd R "--version | head -n 1"
check_cmd git "--version"
echo -e "\nMissing dependencies: $i\n\n\n"

echo "#####################"
echo "#Checking R packages#"
echo "#####################"
Rscript -e \
'needed.packages <- c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI", "BiocManager", "Biostrings")
avail.packages <- installed.packages()
i=0
for( pkg in needed.packages ){ 
    idx <- match(pkg, avail.packages[,"Package"])
    if( ! is.na(idx) ){
        cat("  ",pkg, avail.packages[idx,"Version"], "\n") 
    }else{
        cat("  ", pkg, "NOT FOUND", "\n")
        i=i+1
    }
}
cat("\nMissing R packages: ", i, "\n\n")
'

echo "##############################"
echo "#Checking basic functionality#"
echo "##############################"

Rscript -e \
'
library(sybil)
data(Ec_core)
sybil::optimizeProb(Ec_core)
'
