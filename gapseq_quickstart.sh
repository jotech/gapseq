#!/bin/bash

fasta=$1
model="${fasta%.*}"
media=$2

path=$(readlink -f "$0")
dir=$(dirname "$path")

$dir/./gapseq.sh -b 200 -p all $1

$dir/./transporter.sh $1

Rscript $dir/src/generate_GSdraft.R -r $model-all-Reactions.tbl -t "$model-Transporter.lst" -c $fasta -b 200

Rscript $dir/gf.suite.R -m $model.RDS -n $media -c $model-rxnWeights.RDS -b 50
