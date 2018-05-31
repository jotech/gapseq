#!/bin/bash

curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
echo $dir
metaPwy=$dir/../dat/meta_pwy.tbl
identity=0.9 # clustered uniprot database (0.5 or 0.9)
taxonomy=Bacteria
overwrite=false

# default core metabolism
pwyKey="Amino-Acid-Biosynthesis|Nucleotide-Biosynthesis|Cofactor-Biosynthesis|Carbohydrates-Degradation|CARBO-BIOSYNTHESIS|Polyamine-Biosynthesis|Fatty-acid-biosynthesis|Energy-Metabolism|Terpenoid-Biosynthesis|Chorismate-Biosynthesis"

usage()
{
    echo "Usage"
    echo "$0 -p keyword / -e ec [-t taxonomy]"
    echo "  -p keywords such as pathways or susbstem (default: core metabolism; 'Pathways' for all)"
    echo "  -e search by ec numbers (comma separated)"
    echo "  -t taxonomic range (default: $taxonomy)"
    echo "  -o Should existing files be overwritten (default: $overwrite)"
exit 1
}

while getopts "h?p:e:t:" opt; do
    case "$opt" in
    h|\?)
        usage
        exit 0
        ;;
    p)
        pwyKey=$OPTARG
        ;;
    e)
        ecnumber=$OPTARG
        ;;
    t)
        taxonomy=$OPTARG
        ;;
    esac
done
shift $((OPTIND-1))
[ "$1" = "--" ] && shift

# path for saving sequences
numeric_old=$LC_NUMERIC
LC_NUMERIC="us_US.UTF-8" # must be set in order to get printf working with float numbers
seqpath=$dir/../dat/seq/$taxonomy/unipac$(printf %.0f $(echo "$identity * 100" | bc -l))
LC_NUMERIC=$numeric_old
mkdir -p $seqpath
cd $seqpath

if [ -n "$ecnumber" ]; then
    # create dummpy pwy template for given ec number
    pwyKey=$ecnumber
    pwyDB=$(echo -e "dummy\t$ecnumber\t\t\t\t$ecnumber\t$ecnumber")
else
    # get entries for pathways from database
    pwyDB=$(cat $metaPwy | grep -wEi $pwyKey)
    [ -z "$ecnumber" ] && [ -z "$pwyDB" ] && { echo "No pathways found for key $pwyKey"; exit 1; }
fi

ecs=$(echo "$pwyDB" | awk -F "\t" '{print $7}')
uniq_ecs=$(echo $ecs | tr ',' '\n' | tr ' ' '\n' | sort | uniq)
ecs_max=$(echo "$uniq_ecs" | wc -w)
for ec in $uniq_ecs
do
    i=$((i+1))
    echo -en "\r$i/$ecs_max"
    re="([0-9]+.[0-9]+.[0-9]+.[0-9]+)"
    test=$(if [[ $ec =~ $re ]]; then echo ${BASH_REMATCH[1]}; fi) # check if not trunked ec number (=> too many hits)
    if [ -n "$test" ]; then # check if valid ec
        if [ -f "$ec.fasta" ] &&  [ "$overwrite" = false ]; then # do not update existing files
            continue
        fi
        echo -en " ... Downloading $ec \t\t"
        if [ ! -f "$ec.fasta" ]; then # fasta doesn't exist?
            url="https://www.uniprot.org/uniref/?query=uniprot%3A(ec%3A$ec%20taxonomy%3Abacteria%20AND%20reviewed%3Ayes)%20identity%3A$identity&columns=id%2Creviewed%2Cname%2Ccount%2Cmembers%2Corganisms%2Clength%2Cidentity&format=fasta"
            wget -q $url -O $ec.fasta
        fi
        if [ ! -s "$ec.fasta" ]; then # fasta is empty?
            url="https://www.uniprot.org/uniref/?query=uniprot%3A(ec%3A$ec%20taxonomy%3Abacteria)%20identity%3A$identity&columns=id%2Creviewed%2Cname%2Ccount%2Cmembers%2Corganisms%2Clength%2Cidentity&format=fasta"
            wget -q $url -O $ec.fasta
        fi
    fi
done
