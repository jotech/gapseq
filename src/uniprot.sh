#!/bin/bash

curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
echo $dir
metaPwy=$dir/../dat/meta_pwy.tbl
identity=0.9 # clustered uniprot database (0.5 or 0.9)
identity_unrev=0.5
taxonomy=Bacteria
overwrite=false
get_unrev=false

# default core metabolism
pwyKey="Amino-Acid-Biosynthesis|Nucleotide-Biosynthesis|Cofactor-Biosynthesis|Carbohydrates-Degradation|CARBO-BIOSYNTHESIS|Polyamine-Biosynthesis|Fatty-acid-biosynthesis|Energy-Metabolism|Terpenoid-Biosynthesis|Chorismate-Biosynthesis"

usage()
{
    echo "Usage"
    echo "$0 -p keyword / -e ec [-t taxonomy]"
    echo "  -p keywords such as pathways or susbstem (default: core metabolism; 'Pathways' for all)"
    echo "  -e search by ec numbers (comma separated)"
    echo "  -r search by reaction names (colon separated)"
    echo "  -t taxonomic range (default: $taxonomy)"
    echo "  -o Should existing files be overwritten (default: $overwrite)"
    echo "  -i identity of clustered uniprot database (0.5 or 0.9; default: $identity)"
    echo "  -u get unreviewed instead of reviewed sequences (default $get_unrev)"
    echo "  -g search by gene name"
    echo "  -d search by database reference (e.g. xref)"
exit 1
}

while getopts "h?p:e:t:oi:r:ug:d:" opt; do
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
    o)
        overwrite=true
        ;;
    i)
        identity=$OPTARG
        ;;
    r)
        reaNames=$OPTARG
        ;;
    u)
        get_unrev=true
        ;;
    g)
        geneName=$OPTARG
        ;;
    d)
        dbref=$OPTARG
        ;;
    esac
done
shift $((OPTIND-1))
[ "$1" = "--" ] && shift

# after parsing arguments, only fasta file shoud be there
[ "$#" -ne 0 ] && { usage; }

# consider special case of eg. prokaryotes which are not supported by uniprot taxonomy
if [ "$taxonomy" = "Prokaryota" ]; then
    folder=Prokaryota
    taxonomy="Bacteria,Archaea"
else
    folder=$taxonomy
fi

# path for saving sequences
numeric_old=$LC_NUMERIC
LC_NUMERIC="en_US.UTF-8" # must be set in order to get printf working with float numbers
#seqpath=$dir/../dat/seq/$taxonomy/unipac$(printf %.0f $(echo "$identity * 100" | bc -l))
if [ "$get_unrev" = false ]; then
    seqpath=$dir/../dat/seq/$folder/rev
else
    seqpath=$dir/../dat/seq/$folder/unrev
fi
LC_NUMERIC=$numeric_old
mkdir -p $seqpath
cd $seqpath

[[ ! -f ../updating.log && "$overwrite" = true ]] && echo "date sequences change file status" > ../updating.log

if [ -n "$ecnumber" ]; then
    # create dummpy pwy template for given ec number
    pwyKey=$ecnumber
    pwyDB=$(echo -e "dummy\t$ecnumber\t\t\t\t$ecnumber\t$ecnumber")
elif [ -n "$reaNames" ]; then
    pwyKey=$reaNames
    pwyDB=$(echo -e "dummy\t$reaNames\t\t\t\t$reaNames\t$reaNames")
else
    # get entries for pathways from database
    pwyDB=$(cat $metaPwy | grep -wEi $pwyKey)
    [ -z "$ecnumber" ] && [ -z "$pwyDB" ] && { echo "No pathways found for key $pwyKey"; exit 1; }
fi


if [ -n "$ecnumber" ]; then
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
            [[ -f "$ec.fasta" && "$overwrite" = false ]] && continue # do not update existing files
            [[ ! -f "$ec.fasta" ]] && touch $ec.fasta
            oldmd5=$(md5sum $ec.fasta | cut -d " " -f1)
            oldcount=$(cat $ec.fasta | grep ">" | wc -l)
            rm -f $ec.fasta
            echo -en " ... Downloading $ec \t\t"
            if [ "$get_unrev" = false ]; then
                #swissprot
                Rscript $dir/uniprot_query.R ec $ec $ec.fasta $taxonomy true 0.9 && stat=ok || stat=err
            else
                # unreviewed
                Rscript $dir/uniprot_query.R ec $ec $ec.fasta $taxonomy false 0.5 && stat=ok || stat=err
            fi
            newmd5=$(md5sum $ec.fasta | cut -d " " -f1)
            newcount=$(cat $ec.fasta | grep ">" | wc -l)
            [[ "$oldmd5" != "$newmd5" ]] && echo "`date -u +"%Y-%m-%d_%H:%M:%SUTC"` `echo $newcount` `echo $newcount-$oldcount | bc` $seqpath/$ec.fasta $stat" >> ../updating.log
        fi
    done
fi


if [ -n "$reaNames" ]; then
    reas=$(echo "$pwyDB" | awk -F "\t" '{print $7}')
    uniq_reas=$(echo $reas | tr ';' '\n' | sort | uniq)
    reas_max=$(echo "$uniq_reas" | wc -l)
    for rea in "$uniq_reas"
    do
        i=$((i+1))
        reaNameHash=$(echo -n "$rea" | md5sum | awk '{print $1}')
        echo -en "\r$i/$reas_max"
        [[ -f "$reaNameHash.fasta" && "$overwrite" = false ]] && continue # do not update existing files
        [[ ! -f "$reaNameHash.fasta" ]] && touch $reaNameHash.fasta
        oldmd5=$(md5sum $reaNameHash.fasta | cut -d " " -f1)
        oldcount=$(cat $reaNameHash.fasta | grep ">" | wc -l)
        rm -f $reaNameHash.fasta
        rea="${rea//\"/}"
        echo -en " ... Downloading $rea\t\t"
        if [ "$get_unrev" = false ]; then
            Rscript $dir/uniprot_query.R protein_name "$rea" $reaNameHash.fasta $taxonomy true 0.9 && stat=ok || stat=err
        else
            Rscript $dir/uniprot_query.R protein_name "$rea" $reaNameHash.fasta $taxonomy false 0.5 && stat=ok || stat=err
        fi
        newmd5=$(md5sum $reaNameHash.fasta | cut -d " " -f1)
        newcount=$(cat $reaNameHash.fasta | grep ">" | wc -l)
        [[ "$oldmd5" != "$newmd5" ]] && echo "`date -u +"%Y-%m-%d_%H:%M:%SUTC"` `echo $newcount` `echo $newcount-$oldcount | bc` $seqpath/$reaNameHash.fasta $stat" >> ../updating.log
    done
fi


if [ -n "$geneName" ]; then
    if [ ! -f "$geneName.fasta" ] || [ "$overwrite" = true ]; then # do not update existing files
        [[ ! -f "$geneName.fasta" ]] && touch "$geneName.fasta"
        oldmd5=$(md5sum $geneName.fasta | cut -d " " -f1)
        oldcount=$(cat $geneName.fasta | grep ">" | wc -l)
        rm -f $geneName.fasta
        echo -en " ... Downloading $geneName\t\t"
        if [ "$get_unrev" = false ]; then
            Rscript $dir/uniprot_query.R gene "$geneName" $reaNameHash.fasta $taxonomy true 0.9 && stat=ok || stat=err
        else
            Rscript $dir/uniprot_query.R gene "$geneName" $reaNameHash.fasta $taxonomy false 0.5 && stat=ok || stat=err
        fi
        newmd5=$(md5sum $geneName.fasta | cut -d " " -f1)
        newcount=$(cat $geneName.fasta | grep ">" | wc -l)
        [[ "$oldmd5" != "$newmd5" ]] && echo "`date -u +"%Y-%m-%d_%H:%M:%SUTC"` `echo $newcount` `echo $newcount-$oldcount | bc` $seqpath/$geneName.fasta $stat" >> ../updating.log
    fi
fi

# handling different because several DB links are collected as one reaction sequence file by gapseq_find
if [ -n "$dbref" ]; then
    mkdir -p $seqpath/../rxn
    cd $seqpath/../rxn # different folder necessary because database references are direct sequence links
    if [ ! -f "$dbref.fasta" ] || [ ! "$overwrite" = false ]; then # do not update existing files
        [[ ! -f "$dbref.fasta" ]] && touch "$dbref.fasta"
        oldmd5=$(md5sum $dbref.fasta | cut -d " " -f1)
        oldcount=$(cat $dbref.fasta | grep ">" | wc -l)
        rm -f $dbref.fasta
        echo -en " ... Downloading $dbref\t\t"
        Rscript $dir/uniprot_query.R id $dbref $dbref.fasta $taxonomy && stat=ok || stat=err
        newmd5=$(md5sum $dbref.fasta | cut -d " " -f1)
        newcount=$(cat $dbref.fasta | grep ">" | wc -l)
        [[ "$oldmd5" != "$newmd5" ]] && echo "`date -u +"%Y-%m-%d_%H:%M:%SUTC"` `echo $newcount` `echo $newcount-$oldcount | bc` $seqpath/$dbref.fasta $stat" >> ../updating.log
     fi
fi
