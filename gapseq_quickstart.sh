#!/bin/bash

fasta=$1
if [[ $fasta == *.gz ]]; then # in case fasta is in a archive
    tmp_fasta=$(basename "${fasta}" .gz)
    gunzip -c $fasta > $tmp_fasta
    fasta=$tmp_fasta
fi
model="${fasta%.*}"
media=$2

path=$(readlink -f "$0")
dir=$(dirname "$path")

$dir/./gapseq.sh -b 200 -p all $fasta

$dir/./transporter.sh $fasta

Rscript $dir/src/generate_GSdraft.R -r $model-all-Reactions.tbl -t "$model-Transporter.lst" -c $fasta -b 200

Rscript $dir/gf.suite.R -m $model.RDS -n $media -c $model-rxnWeights.RDS -b 50
