#!/bin/bash

curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
metaRea=$dir/../dat/meta_rea.tbl
reaDB1=$dir/../dat/vmh_reactions.tsv
reaDB2=$dir/../dat/bigg_reactions.tbl
reaDB3=$dir/../dat/seed_reactions_corrected.tsv
reaDB4=$dir/../dat/mnxref_seed.tsv
reaDB5=$dir/../dat/mnxref_seed-other.tsv
reaDB6=$dir/../dat/mnxref_bigg-other.tsv
brenda=$dir/../dat/brenda_ec_edited.csv
seedEC=$dir/../dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv
seedEnzymesNames=$dir/../dat/seed_Enzyme_Name_Reactions_Aliases.tsv
altecdb=$dir/../dat/altec.csv


rea=$1
reaName=$2
ec=$3
database=$4
EC_test=$5

[ $# -eq 0 ] && { echo "Usage: $0 id name ec database ec_test"; exit 1; }

#echo rea=$1 reaName=$2 ec=$3 database=$4 ec_test=$5

getDBhit(){
    [[ -n "$rea" ]] && kegg=$(grep -wFe "|$rea" $metaRea | awk -F "\t" {'print $5'})

    # 1) search in reaction db by EC
    if [[ -n "$EC_test" ]]; then
        if [ "$database" == "vmh" ]; then
            dbhit=$(grep -wF $ec $reaDB1 | awk -F ',' '{print $1}')
        elif [ "$database" == "seed" ]; then
            dbhit=$(cat $seedEC | cut -f1,3 | grep -wF $ec | cut -f1 | tr '|' ' ')
            #dbhit="$dbhit $(grep -wF $ec $reaDB4 | awk -F '\t' '{print $4}' | tr '\n' ' ')" # mnxref considers also obsolete seed hits 
        fi
    fi

    # 2) search in reaction db by kegg identifier 
    if [ "$database" == "vmh" ]; then
        [ -n "$kegg" ]  && dbhit="$dbhit $(grep -wE "$(echo $kegg |tr ' ' '|')" $reaDB1 | awk -F ',' '{print $1}')"
    elif [ "$database" == "seed" ]; then
        [ -n "$kegg" ] && dbhit="$dbhit $(grep -wE "$(echo $kegg | tr ' ' '|')" $reaDB3 | awk -F '\t' '$18 == "OK" {print $1}' )" # only consider reactions which are OK
        #[ -n "$kegg" ] && dbhit="$dbhit $(grep -wE "$(echo $kegg | tr ' ' '|')" $reaDB4 | awk -F '\t' '{print $4}' | tr '\n' ' ')"
        [ -n "$kegg" ] && dbhit="$dbhit $(grep -wE "$(echo $kegg | tr ' ' '|')" $reaDB5 | awk -F '\t' '{print $2}' | tr '\n' ' ')" # use single database (reaDB4,5 are similiar)
    fi

    # 3) search in reaction db by alternative EC
    if [[ -n "$EC_test" ]]; then
        altec=$(grep -wF $ec $altecdb | grep -P "([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)" -o | grep -v $ec)
        if [ -z "$altec" ]; then  
            brendaec=$(grep -wF $ec $brenda | grep -P "([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)" -o | grep -v $ec)
            [[ `echo "$brendaec" | wc -l` -eq 1 ]] && altec=brendaec # only take unique hits (too many multiple transferred ECs)
        fi

        if [ "$database" == "vmh" ]; then
            [ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo ${altec//./\\.} | tr ' ' '|')" $reaDB1 | awk -F ',' '{print $1}')" # take care of multiple EC numbers
        elif [ "$database" == "seed" ]; then
            [ -n "$altec" ] && dbhit="$dbhit $(cat $seedEC | cut -f1,3 | grep -wE "$(echo ${altec//./\\.} | tr ' ' '|')" | cut -f1 | tr '|' ' ')" # take care of multiple EC numbers
            #[ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo ${altec//./\\.} | tr ' ' '|')" $reaDB4 | awk -F '\t' '{print $4}')" # mnxref considers also obsolete seed hits 
        fi
    fi

    # 4) search in bigg db by metacyc id (does only make sense for vmh/bigg namespace)
    if [ "$database" == "vmh" ] && [ -n "$rea" ]; then
        dbhit="$dbhit $(grep -wFe "$rea" $reaDB2 | awk '{print $1}')"
    fi

    # 5) match reaction using mnxref namespace
    if [ "$database" == "seed" ] && [ -n "$rea" ]; then
        dbhit="$dbhit $(grep -wFe "|$rea" $reaDB5 | awk '{print $2}')"
    elif [ "$database" == "vmh" ] && [ -n "$rea" ]; then
        dbhit="$dbhit $(grep -wFe "|$rea" $reaDB6 | awk '{print $2}')"
    fi

    # 6) match reaction using custom enzyme-name - seedID mapping
    if [ "$database" == "seed" ] && [ -n "$reaName" ]; then
        dbhit="$dbhit $(grep -wFe "$reaName" $seedEnzymesNames | awk -F '\t' ' {print $1}')"
    fi

    [ "$dbhit" == " " ] && dbhit=""

}

if [[ -z "$EC_test" ]]; then
    re="([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)"
    EC_test=$(if [[ $ec =~ $re ]]; then echo ${BASH_REMATCH[1]}; fi) # check if not trunked ec number (=> too many hits)
fi
[[ -z "$database" ]] && database="seed"

getDBhit
dbhit="$(echo $dbhit | tr ' ' '\n' | sort | uniq | tr '\n' ' ')"

echo $dbhit
