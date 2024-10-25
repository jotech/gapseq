#!/bin/bash

usage()
{
    echo "$0 [nr] [fast]"
    echo "   nr    certain line number of organism in model_flux.tbl (optional, otherwise all organism are considered)"
    echo "   fast  use old files and only redo gap gilling step (optional, by default all steps will be performed)"
    exit 0
}
[[ "$#" -gt 2 ]] && usage

curdir=$(pwd)                                                                   
path=$(readlink -f "$0")
dir=$(dirname "$path")
db=$dir/model_flux.tbl
db_rows=$(cat $db | wc -l)
version=$($dir/../gapseq -v | sed 's/gapseq version: //')
log_file="$dir/model_flux.log"
fast=false

echo $0 $dir
mkdir -p $dir/genomes
mkdir -p $dir/out
cd $(mktemp -d)

if [[ -n "$1" ]] && [[ ! "$1" == "fast" ]]; then
    [[ ( $1 -gt $db_rows ) || ( $1 -le 1 ) ]] && { echo "Line $1 does not exist in database (range: 2-$db_rows)"; exit 1; }
    range=$1
else
    range=$(seq $db_rows) 
fi
if [[ "$1" == "fast" ]] || [[ "$2" == "fast" ]]; then
    fast=true
    echo Fast mode, only gap filling is done.
fi

#download genomes
wget -q --spider https://www.ncbi.nlm.nih.gov
if [[ $? -eq 0 ]]; then
    needed=$(mktemp)
    avail=$(mktemp)
    cat $dir/model_flux.tbl | cut -f2 | tail -n +2 | sort > $needed
    keywords=$(cat $needed | tr '\n' '|' | rev | cut -c2- | rev)
    ls $dir/genomes | grep -Eo $keywords |sort > $avail
    download=$(comm -23 $needed $avail | tr '\n' ',' | rev | cut -c2- | rev)
    if [[ -n "$download" ]]; then
        echo Downloading $download
        ncbi-genome-download -A $download --flat-output -F protein-fasta -o $dir/genomes/ bacteria
    fi
fi

for i in $range
do
    [[ $i -eq 1 ]] && continue
    line=$(sed -n ${i}p $db)
    name=$(echo "$line" | cut -f1)
    genome=$(echo "$line" | cut -f2)
    genome_file=$(ls $dir/genomes/${genome}*.faa.gz)
    media=$(echo "$line" | cut -f3 | tr ',' '\n')
    base=$(basename "$genome_file")
    rxn=$(echo "$line" | cut -f4)
    rxn_dir=$(echo "$line" | cut -f5)
    description=$(echo "$line" | cut -f6)
    parm_find=$(echo "$line" | cut -f7)
    parm_transport=$(echo "$line" | cut -f8)
    parm_draft=$(echo "$line" | cut -f9)
    parm_fill=$(echo "$line" | cut -f10)
    
    id_tmp=${base%.*}
    id=${id_tmp%.*}
    echo $i $name $id $media

    if [[ "$fast" = false ]]; then
        $dir/../gapseq find $parm_find -v 0 -p all "$genome_file"
        mv "${id}-all-Reactions.tbl" "${id}-all-Pathways.tbl" $dir/out/
        $dir/../gapseq find-transport $parm_transport "$genome_file"
        mv "${id}-Transporter.tbl" $dir/out/
        $dir/../gapseq draft $parm_draft -r $dir/out/"$id-all-Reactions.tbl" -t $dir/out/"$id-Transporter.tbl" -c "$genome_file" -p "$id-all-Pathways.tbl"
        mv "${id}-draft.RDS" "${id}-rxnWeights.RDS" "${id}-rxnXgenes.RDS" $dir/out/
    fi
    #ls "$dir/out/${id}-all-Reactions.tbl" "$dir/out/${id}-all-Pathways.tbl" "$dir/out/${id}-Transporter.tbl" "$dir/out/${id}-draft.RDS" "$dir/out/${id}-rxnWeights.RDS" "$dir/out/${id}-rxnXgenes.RDS"

    for medium in $media
    do
        medium_file=$(ls $dir/../dat/media/${medium}.csv)
        $dir/../gapseq fill $parm_fill -m $dir/out/"${id}-draft.RDS" -n "$medium_file" -c $dir/out/"${id}-rxnWeights.RDS" -g $dir/out/"${id}-rxnXgenes.RDS"
        mv ${id}.RDS $dir/out/${id}_${medium}.RDS
        model_file=$(ls $dir/out/${id}_${medium}.RDS)

        Rscript -e '
        suppressMessages(library(cobrar))
        suppressMessages(library(stringr))

        args = commandArgs(trailingOnly=TRUE)
        mod.file <- args[1]
        rxn <- unlist(str_split(args[2], ","))
        dir <- unlist(str_split(args[3], ","))
        description <- args[4]
        org.name <- args[5]
        version  <- args[6]
        medium   <- args[7]
        log.file <- args[8]

        mod <- readRDS(mod.file)
        sol <- fba(mod)
        fva <- fva(mod, react=rxn)
        fva_max <- fva$max.flux
        fva_min <- fva$min.flux

        cat("Checking phenotype:", description, "on medium:", medium, "\n")
        rxn.all <- 0
        for(i in seq_along(rxn)){
            rxn.id  <- rxn[i]
            rxn.idx <- match(rxn.id, mod@react_id)
            if( is.na(rxn.idx) ) warning(paste("Reaction", rxn[i], "not found in model", mod.file))
            rxn.dir <- dir[i]
            rxn.sol <- sol@fluxes[rxn.idx]
            rxn.max <- fva_max[i]
            rxn.min <- fva_min[i]
            if( rxn.dir==">" ){
                if( rxn.sol>0 | rxn.max > 0) rxn.cor <- TRUE else rxn.cor <- FALSE
                }else if( rxn.dir=="<" ){
                if( rxn.sol<0 | rxn.min < 0) rxn.cor <- TRUE else rxn.cor <- FALSE
                } else {
                Stop(paste("Direction not defined for reaction:", rxn.id, rxn.dir))
            }
            rxn.all <- rxn.all + rxn.cor
            cat("  ", rxn.id, "predicted:", round(rxn.sol,1), "(",round(rxn.min,1),",",round(rxn.max,1),")", "expected:", rxn.dir, "Status:", rxn.cor, "\n")
        }
        all.cor <- round(rxn.all*100/length(rxn),1)
        cat("Phenotyp correctness:", all.cor, "%\n\n")
        line=paste(Sys.Date(), version, basename(mod.file), org.name, description, medium, all.cor, sep="\t")
        write(line, file=log.file, append=TRUE)
        ' $model_file "$rxn" "$rxn_dir" "$description" "$name" "$version" "$medium" "$log_file"
    done
done

cd $curdir
exit 0
