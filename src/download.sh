#!/bin/bash


curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
seqpath=$dir/../dat/seq
echo Target directory: $seqpath

echo Downloading uniprot swisspot database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

gunzip uniprot_sprot.fasta.gz
mv uniprot_sprot.fasta $seqpath

echo Building blast database 
makeblastdb -in $seqpath/uniprot_sprot.fasta -dbtype prot -out $seqpath/uniprot_sprot
