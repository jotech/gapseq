#!/bin/bash

# This small script checks a fasta file (uncompressed) if it is a protein fasta ("prot") or
# nucleotide fasta ("nucl").

input_mode=""

n_char=`cat $1 | grep -v "^>" | awk '{for(i=1;i<=NF;i++)if(!a[$i]++)print $i}' FS="" | wc -l`
if [ $n_char -ge 15 ]; then
    input_mode="prot"
else
    input_mode="nucl"
fi

echo $input_mode

