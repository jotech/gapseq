#!/bin/bash

# Check if pyrodigal is installed
if ! command -v pyrodigal >/dev/null 2>&1; then
  echo "pyrodigal is NOT installed, but required to translate nucleotide genome to protein amino acid sequences."
  exit 1
fi


curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")

codes=(4 11 25)
outpfx=""
infna=""
n_threads=1
verbose=1

usage()
{
    echo "Usage"
    echo "$0 -i fasta with genomic sequence(s)"
    echo "  -o Prefix for file exports"
    echo "  -c Genetic code (i.e., BLAST translation table Id)"
    echo "  -K Number of threads to use"
exit 1
}

while getopts "h?p:i:o:c:K:v:" opt; do
    case "$opt" in
    h|\?)
        usage
        exit 0
        ;;
    i)
        infna=$OPTARG
        ;;
    o)
        outpfx=$OPTARG
        ;;
    c)
        code_tmp=$OPTARG
        if [ $code_tmp != "auto" ]; then
            codes=$code_tmp
        fi
        ;;
    K)
        n_threads=$OPTARG
        ;;
    esac
done
shift $((OPTIND-1))
[ "$1" = "--" ] && shift

for code in "${codes[@]}"; do
    pyrodigal -i $infna -a ${outpfx}_${code}.faa -g $code -j $n_threads -o ${outpfx}_${code}.gff
done

if [ "${#codes[@]}" -gt 1 ]; then
    Rscript $dir/predict_codontable.R $infna
    code=`cat ${outpfx}_code`
else
    code=$codes
fi

mv "${outpfx}_${code}.faa" "${outpfx}.faa"
mv "${outpfx}_${code}.gff" "${outpfx}.gff"

