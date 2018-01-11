#!/bin/bash

# TODO: dbhit: handle when $ec has more then one number
# TODO: handle incomplete/unspecific ecs from metacyc (e.g. get ec from kegg, update maually or get genes from metacyc)
# TODO: advanced handling of parameters/arguments (argv), e.g. output file)


# paths and variables
sbml=$(readlink -f $3)
fasta=$(readlink -f $2)
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
seqpath=$dir/dat/seq/
bitcutoff=50 # cutoff blast: min bit score
metaPwy=$dir/dat/meta_pwy.tbl
metaRea=$dir/dat/meta_rea.tbl
reaDB1=$dir/dat/vmh_reactions.csv
reaDB2=$dir/dat/bigg_reactions.tbl
brenda=$dir/dat/brenda_ec.csv


# set databases
[ "$1" == "all" ]       && pwyKey=Pathways
[ "$1" == "amino" ]     && pwyKey=Amino-Acid-Biosynthesis
[ "$1" == "nucl" ]      && pwyKey=Nucleotide-Biosynthesis
[ "$1" == "cofactor" ]  && pwyKey=Cofactor-Biosynthesis
[ "$1" == "carbo" ]     && pwyKey=Carbohydrates-Degradation
[ "$1" == "polyamine" ] && pwyKey=Polyamine-Biosynthesis
[ "$1" == "fatty" ]     && pwyKey=Fatty-acid-biosynthesis
[ "$1" == "core" ]      && pwyKey="Amino-Acid-Biosynthesis|Nucleotide-Biosynthesis|Cofactor-Biosynthesis|Carbohydrates-Degradation|Polyamine-Biosynthesis|Fatty-acid-biosynthesis"

[ -z $pwyKey ] && pwyKey=$1

# USAGE
#[ $# -ne 2 ] && { echo "Usage: $0 file.fasta model.sbml"; exit 1; }
( [ $# -ne 3 ] ) && { echo "Usage: $0 database (amino,nucl,cofactor,carbo,polyamine) file.fasta model.sbml"; exit 1; }

pwyDB=$(cat $metaPwy | grep -wE $pwyKey)
[ -z "$pwyDB" ] && { echo "No pathways found for key $pwyKey"; exit 1; }


# tmp working directory
cd $(mktemp -d)
#cd /tmp/tmp.VMnit0ThVM 


# try to read sbml file and save as r object
$dir/src/sbml_read.R $sbml


# create blast database
makeblastdb -in $fasta -dbtype nucl -out orgdb >/dev/null


cand=""     #list of candidate reactions to be added
bestCand="" # list of candidates from (almost) complete pathways
bestPwy=""  # list of found pathways

echo Checking for reaction from: $1 $pwyKey
for i in `seq 1 $(echo "$pwyDB" | wc -l)`
do
    pwyCand="" # candidate reaction of current pathway
    count=0
    countex=0
    countdb=0
    vague=0
    line=$(echo "$pwyDB" | awk -v i=$i 'NR==i')
    pwy=$(echo "$line" | awk -F "\t" '{print $1}')
    name=$(echo "$line" | awk -F "\t" '{print $2}')
    ecs=$(echo "$line" | awk -F "\t" '{print $7}')
    reaids=$(echo "$line" | awk -F "\t" '{print $6}')
    echo -e '\n'$i. Checking for pathway $pwy $name
    for j in `seq 1 $(echo $ecs | tr "," "\n"i | wc -l)`
    #for ec in $(echo $ecs | tr "," "\n")
    do 
        ec=$(echo $ecs | awk -v j=$j -F ',' '{print $j}')
        rea=$(echo $reaids | awk -v j=$j -F ',' '{print $j}')
        re="([0-9]+.[0-9]+.[0-9]+.[0-9]+)"
        test=$(if [[ $ec =~ $re ]]; then echo ${BASH_REMATCH[1]}; fi) # check if not trunked ec number (=> too many hits)
        if [ -n "$test" ]; then
            ((count++))
            query=$seqpath$ec.fasta
            if [ -s $query ]; then
                out=$pwy-$ec.blast
                tblastn -db orgdb -query $query -outfmt '6 qseqid sseqid pident evalue bitscore stitle' >$out 
                if [ -s $out ]; then
                    bhit=$(cat $out | awk -v bitcutoff=$bitcutoff '$5>=bitcutoff {print $1}' | wc -l)
                    if [ $bhit -gt 0 ]; then
                        echo -e '\t'Blast hit: $rea $ec
                        
                        
                        kegg=$(grep -wF $rea $metaRea | awk -F "\t" {'print $5'})
                        # search in vmh db by EC
                        altec=$(grep $ec $brenda | grep -P "([0-9]+.[0-9]+.[0-9]+.[0-9]+)" -o | grep -v $ec)
                        # search in vmh db by kegg identifier 
                        dbhit=$(grep -wF $ec $reaDB1 | awk -F ',' '{print $1}')
                        [ -n "$kegg" ]  && [ -z "$dbhit" ]&& dbhit=$(grep -wF $kegg $reaDB1 | awk -F ',' '{print $1}')
                        # search in vmh db by alternative EC
                        [ -n "$altec" ] && [ -z "$dbhit" ]  && dbhit=$(grep -wE "$(echo $altec | tr ' ' '|')" $reaDB1 | awk -F ',' '{print $1}') # take care of multiple EC numbers
                        # search in bigg db by metacyc id
                        [ -z "$dbhit" ]  && dbhit=$(grep -wF $rea $reaDB2 | awk '{print $1}')
                        if [ -n "$dbhit" ]; then
                            echo -e '\t\t'Candidate reaction for import: $dbhit
                            pwyCand="$pwyCand$dbhit " # remember candidate reaction
                            ((countdb++))
                        else
                            echo -e '\t\t'No candidate reaction found for import: $rea $ec
                        fi
                        ((countex++))
                    else
                        echo -e '\t'No significant blast hits found for reaction: $rea $ec
                    fi
                else
                    echo -e '\t'No blast hits found for reaction: $rea $ec
                fi
            else
                echo -e '\t'No sequence data found for $rea $ec ..skipping..
               ((vague++))
            fi
        else
            #TODO: if unspecific should the reaction still be added?
            echo -e '\t'EC number too unspecific: $rea $ec ..skipping..
            ((vague++))
        fi
    done
    echo -e Pathway completness: $countex/$count with $vague reactions of unclear state
    echo -e Hits with candidate reactions in database: $countdb/$count
    cand="$cand$pwyCand "
    if [ $count -ne 0 ] && [ $(echo "scale=2; $countex/$count>=0.9" | bc) -ne 0 ]; then
        bestCand="$bestCand$pwyCand " # save candidates from almost complete (>=90%) pathways
        bestPwy="$bestPwy$name\n"
    fi
done

cand=$(echo $cand | tr ' ' '\n' | sort | uniq | tr '\n' ' ') # remove duplicates
echo -e '\n'Total candidate reactions:
echo $cand

echo -e '\n'Pathways found:
echo -e $bestPwy

bestCand=$(echo $bestCand | tr ' ' '\n' | sort | uniq | tr '\n' ' ') # remove duplicates
echo -e '\n'Candidate reactions from complete pathways:
echo -e $bestCand

# add reactions and write new sbml model
echo $bestCand > newReactions.lst
echo ""
$dir/src/sbml_write.R newReactions.lst $dir
modelold=$(basename $sbml)
modelnew="${modelold%.*}G.xml"
cp modelnew.xml $curdir/$modelnew
