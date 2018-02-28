#!/bin/bash

# need barrnap https://github.com/tseemann/barrnap
# bedtools https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html

echo Takes a fasta file, extracts 16S rRNA, and does blast search against agora 16S db
# show usage
[ $# -eq 0 ] && { echo "Usage: $0 file.fasta"; exit 1; }

path=$(readlink -f "$0")
dir=$(dirname "$path")
input=$1
db=$dir/../dat/agora16S
tmp_gff=`mktemp`
tmp_fasta=`mktemp`

# find 16S
printf "\nTry to find 16S sequences\n\n"
barrnap -q $input | grep 16S > $tmp_gff
cat $tmp_gff
[ ! -s $tmp_gff ] && { echo "No 16S rRNA found"; exit 1; }

# extract 16S
bedtools getfasta -fi $input -bed $tmp_gff -fo $tmp_fasta

# blast against ncbi 16S db
printf "\nBlast against agora 16S database\n\n"
blastn -db $db -query $tmp_fasta -perc_identity 97.5 -outfmt '6 sseqid pident evalue bitscore stitle'

# delete index file created by bedtools
rm $input.fai
