#!/bin/bash

echo $1

curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")

# get path to ldconfig (missing in e.g. debian)
ldconfig2=$(whereis -b ldconfig | cut -f2 -d " ")

# temporary working directory
cd $(mktemp -d)

# gapseq version
echo -e `$dir/../gapseq -v | head -n 1`

# printing operation system
echo $OSTYPE
echo -e `uname -v` "\n\n"

echo "#######################"
echo "#Checking dependencies#"
echo "#######################"
i=0

check_cmd(){
    cmd=$1
    arg=$2
    print=$3
    altname=$4
    [[ -n "$altname" ]] && name=$altname || name=$cmd
    if command -v $cmd &> /dev/null
    then
        version=$(eval $cmd "$arg")
        if [[ -n "$version" ]]; then
            [[ ! "$print" = false ]] && echo $version
        else
            echo $name NOT FOUND
            i=$((i+1))
        fi
    else
        echo $name NOT FOUND
        i=$((i+1))
    fi
}

# file-based library check for mac os
get_lib_ver_macos() { otool -L "$1" 2>/dev/null | sed -n '2p' | sed -E 's/.*current version ([^,)]+).*/\1/'; }

check_lib_macos() {
    libname=$1
    for d in /usr/lib /usr/local/lib /opt/homebrew/lib $CONDA_PREFIX/lib ${LD_LIBRARY_PATH//:/ }; do
        if [[ -d "$d" ]]; then
            libfile=$(find -L "$d" -name "${libname}.dylib" -type f -maxdepth 1 2>/dev/null | head -n 1)
            if [[ -n "$libfile" ]]; then
                version=$(get_lib_ver_macos "$libfile")
                echo "$libname ${version:-found} ($libfile)"
                return 0
            fi
        fi
    done
    echo "$libname NOT FOUND"
    i=$((i+1))
}

if [[ "$OSTYPE" == "darwin"* ]]; then
    check_lib_macos "libsbml"
    check_lib_macos "libglpk"
    check_cmd awk "--version | head -n 1" # one true awk (macos/bsd)
else
    check_cmd $ldconfig2 "-V | head -n 1"
    check_cmd $ldconfig2 "-N -v $(sed 's/:/ /g' <<< $LD_LIBRARY_PATH:$CONDA_PREFIX/lib) 2>/dev/null | grep sbml.so" true libsbml
    check_cmd $ldconfig2 "-N -v $(sed 's/:/ /g' <<< $LD_LIBRARY_PATH:$CONDA_PREFIX/lib) 2>/dev/null | grep glpk.so" true libglpk
    check_cmd awk "-W version 2>/dev/null | head -n 1" # works with both GNU awk and mawk
fi
check_cmd sed "--version | head -n 1"
check_cmd grep "-V | head -n 1"
check_cmd blastp "-version | head -n 1"
check_cmd R "--version | head -n 1"
check_cmd Rscript "--version | head -n 1" false
check_cmd git "--version"
check_cmd hmmsearch "-h | grep \"# HMMER\" | cut -c 3-"
check_cmd bc "-v | head -n 1"
echo -e "\nMissing dependencies: $i\n\n"

echo "#####################"
echo "#Checking R packages#"
echo "#####################"
Rscript -e 'needed.packages <- c("data.table", "stringr", "cobrar", "getopt", "R.utils", "stringi", "BiocManager", "Biostrings", "jsonlite", "httr"); avail.packages <- installed.packages(); i=0; for( pkg in needed.packages ){; idx <- match(pkg, avail.packages[,"Package"]); if( ! is.na(idx) ){; cat(pkg, avail.packages[idx,"Version"], "\n"); }else{; cat(pkg, "NOT FOUND", "\n"); i=i+1; }; }; cat("\nMissing R packages: ", i, "\n\n\n")'
Rscript -e 'if( packageVersion("cobrar") < "0.1.2" ) cat("WRONG cobrar version (>0.1.2 needed)\n\n\n")'

echo "##############################"
echo "#Checking basic functionality#"
echo "##############################"
i=0

# Optimization
Rscript -e 'suppressMessages(library(cobrar)); fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar"); mod <- readSBMLmod(fpath); sol <- fba(mod); sol.stat <- sol@stat %in% c(2,5) && sol@obj > 1e-7; cat("Optimization test:", ifelse( sol.stat, "OK", "FAILED"), "\n"); ' > opt.log
cat opt.log
if  grep -q "OK" opt.log; then
    i=$((i+1))
fi

# Construct full model
Rscript -e 'args = commandArgs(trailingOnly=TRUE); source(paste0(args[1],"/construct_full_model.R")); fmod = construct_full_model(args[1]); cat("Building full model:",ifelse(class(fmod) == "ModelOrg", "OK", "FAILED"), "\n"); ' $dir > opt.log
cat opt.log
if  grep -q "OK" opt.log; then
    i=$((i+1))
fi

# blast
fasta=$dir/../dat/seq/Bacteria/user/1.2.1.87.fasta
tmp_file=$(mktemp)
gunzip -c $dir/../toy/myb71.faa.gz > $tmp_file
makeblastdb -in $tmp_file -dbtype prot -out orgdb >/dev/null
blastp -db orgdb -query $fasta -qcov_hsp_perc 75 -outfmt "6 $blast_format" > blast.out
if [[ -s blast.out ]]; then
    echo "Blast test: OK"
    i=$((i+1))
else
    echo "Blast test: FAILED"
fi

echo -e "\nPassed tests: $i/3\n"

echo "################################"
echo "#Checking optional dependencies#"
echo "################################"

echo "(Alternative sequence alignment tools)"
echo "    `check_cmd diamond "--version "`"
echo "    mmseqs `check_cmd mmseqs "version"`"

echo "(Nucleotide genome to protein genome translation)"
echo "    `check_cmd pyrodigal "--version "`"
echo ""

