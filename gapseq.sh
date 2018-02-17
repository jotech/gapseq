#!/bin/bash

# TODO: limit pwy search to taxonomic scope
# TODO: save dummy seq file for ec without uniprot hit (save nonsense requests)
# TODO: handle incomplete/unspecific ecs from metacyc (e.g. get ec from kegg, update maually or get genes from metacyc)
# TODO: if taxonomic range is not bacteria, then sequence data must be updated!

pathways=""
input_sbml=""
output_sbml=""
database="seed"
verbose=0
taxonomy="Bacteria"
bitcutoff=50 # cutoff blast: min bit score
identcutoff=0   # cutoff blast: min identity
covcutoff=75 # cutoff blast: min coverage
strictCandidates=false
completnessCutoff=66 # consider pathway to be present if other hints (e.g. key enzyme present) are avaiable and pathway completness is at least as high as completnessCutoff (requires strictCandidates=false)

usage()
{
    echo "Usage"
    echo "$0 -p keyword [-d database] [-f model.sbml] [-t taxonomy] file.fasta."
    echo "  -p keywords such as pathways or susbstem (for example amino,nucl,cofactor,carbo,polyamine)"
    echo "  -d database: vmh or seed (default: $database)"
    echo "  -t taxonomic range (default: $taxonomy)"
    echo "  -b bit score cutoff for local alignment (default: $bitcutoff)"
    echo "  -i identity cutoff for local alignment (default: $identcutoff)"
    echo "  -c coverage cutoff for local alignment (default: $covcutoff)"
    echo "  -s strict candidate reaction handling (do _not_ use pathway completness, key kenzymes and operon structure to infere if imcomplete pathway could be still present (default: $strictCandidates)"
exit 1
}
# USAGE
#[ $# -ne 2 ] && { echo "Usage: $0 file.fasta model.sbml"; exit 1; }
#( [ $# -lt 2 ] || [ $# -gt 3 ] ) && { usage; }
#( [ $# -eq 3 ] ) && { withSbml=true; sbml=$(readlink -f $3); }


# paths and variables
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
seqpath=$dir/dat/seq/
metaPwy=$dir/dat/meta_pwy.tbl
metaRea=$dir/dat/meta_rea.tbl
reaDB1=$dir/dat/vmh_reactions.csv
reaDB2=$dir/dat/bigg_reactions.tbl
reaDB3=$dir/dat/seed_reactions.tsv
reaDB4=$dir/dat/mnxref_seed.tsv
brenda=$dir/dat/brenda_ec.csv
seedEC=$dir/dat/seed_Enzyme_Class_Reactions_Aliases_unique.tsv


# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "h?p:d:i:b:c:vs:o:t:s" opt; do
    case "$opt" in
    h|\?)
        usage
        exit 0
        ;;
    p)  
        pathways=$OPTARG
        ;;
    d)  
        database=$OPTARG
        ;;
    v)  
        verbose=1
        ;;
    s)  
        input_sbml=$(readlink -f $OPTARG)
        $dir/src/sbml_read.R $input_sbml
        ;;
    b)
        bitcutoff=$OPTARG
        ;;
    i)
        identcutoff=$OPTARG
        ;;
    c)
        covcutoff=$OPTARG
        ;;
    t)  
        taxonomy=$OPTARG
        ;;
    o)  
        output_sbml=$(readlink -f $OPTARG)
        ;;
    s)
        strictCandidates=true
    esac
done
shift $((OPTIND-1))
[ "$1" = "--" ] && shift

# after parsing arguments, only fasta file shoud be there
[ "$#" -ne 1 ] && { usage; }

# get fasta file
fasta=$(readlink -f "$1")
[ ! -s $fasta ] && { echo Invalid file: $1; exit 0; }
tmpvar=$(basename $fasta)
fastaID="${tmpvar%.*}"

# pathways and fasta file have to be provided
( [ -z "$pathways" ] || [ -z "$fasta" ] ) && { usage; }

# select pathway keys to be used in database search
case $pathways in
    all)
        pwyKey=Pathways
        ;;
    amino)
        pwyKey=Amino-Acid-Biosynthesis
        ;;
    nucl)
        pwyKey=Nucleotide-Biosynthesis
        ;;
    cofactor)
        pwyKey=Cofactor-Biosynthesis
        ;;
    carbo)
        pwyKey=CARBO-BIOSYNTHESIS
        ;;
    carbo-deg)
        pwyKey=Carbohydrates-Degradation
        ;;
    polyamine)
        pwyKey=Polyamine-Biosynthesis
        ;;
    fatty)
        pwyKey=Fatty-acid-biosynthesis
        ;;
    energy)
        pwyKey=Energy-Metabolism
        ;;
    terpenoid)
        pwyKey=Terpenoid-Biosynthesis
        ;;
    degradation)
        pwyKey=Degradation
        ;;
    core)
        pwyKey="Amino-Acid-Biosynthesis|Nucleotide-Biosynthesis|Cofactor-Biosynthesis|Carbohydrates-Degradation|CARBO-BIOSYNTHESIS|Polyamine-Biosynthesis|Fatty-acid-biosynthesis|Energy-Metabolism|Terpenoid-Biosynthesis"
        ;;
    *)
        pwyKey=$pathways
        ;;
esac


pwyDB=$(cat $metaPwy | grep -wE $pwyKey)
[ -z "$pwyDB" ] && { echo "No pathways found for key $pwyKey"; exit 1; }



# function to get database hits for ec number
getDBhit(){
    kegg=$(grep -wF $rea $metaRea | awk -F "\t" {'print $5'})
    altec=$(grep $ec $brenda | grep -P "([0-9]+.[0-9]+.[0-9]+.[0-9]+)" -o | grep -v $ec)
    
    # 1) search in reaction db by EC
    if [ "$database" == "vmh" ]; then
        dbhit=$(grep -wF $ec $reaDB1 | awk -F ',' '{print $1}')
    elif [ "$database" == "seed" ]; then
        dbhit=$(grep -wF $ec $seedEC | awk '{print $1}' | tr '|' '\n')
        dbhit="$dbhit $(grep -wF $ec $reaDB4 | awk -F '\t' '{print $4}')"
    fi

    # 2) search in reaction db by kegg identifier 
    if [ "$database" == "vmh" ]; then
        [ -n "$kegg" ]  && dbhit="$dbhit $(grep -wE "$(echo $kegg |tr ' ' '|')" $reaDB1 | awk -F ',' '{print $1}')"
    elif [ "$database" == "seed" ]; then
        [ -n "$kegg" ] && dbhit="$dbhit $(grep -wE "$(echo $kegg | tr ' ' '|')" $reaDB3 | awk -F '\t' '$18 == "OK" {print $1}' )" # only consider reactions which are OK
        [ -n "$kegg" ] && dbhit="$dbhit $(grep -wE "$(echo $kegg | tr ' ' '|')" $reaDB4 | awk -F '\t' '{print $4}')"
    fi

    # 3) search in reaction db by alternative EC
    if [ "$database" == "vmh" ]; then
        [ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo $altec | tr ' ' '|')" $reaDB1 | awk -F ',' '{print $1}')" # take care of multiple EC numbers
    elif [ "$database" == "seed" ]; then
        [ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo $altec | tr ' ' '|')" $seedEC | awk '{print $1}' | tr '|' '\n')" # take care of multiple EC numbers
        [ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo $altec | tr ' ' '|')" $reaDB4 | awk -F '\t' '{print $4}')" # take care of multiple EC numbers
    fi
    
    # 4) search in bigg db by metacyc id (does only make sense for vmh/bigg namespace)
    if [ "$database" == "vmh" ]; then
        dbhit="$dbhit $(grep -wF $rea $reaDB2 | awk '{print $1}')"
    fi

    [ "$dbhit" == " " ] && dbhit=""
}




# tmp working directory
cd $(mktemp -d)

# create blast database
makeblastdb -in $fasta -dbtype nucl -out orgdb >/dev/null


cand=""     #list of candidate reactions to be added
bestCand="" # list of candidates from (almost) complete pathways
bestPwy=""  # list of found pathways

pwyNr=$(echo "$pwyDB" | wc -l)
echo Checking for pathways and reactions in: $1 $pwyKey
echo Number of pathways to be considered: $pwyNr
for i in `seq 1 $pwyNr`
do
    pwyCand="" # candidate reaction of current pathway
    pwyCandAll="" # all possible reaction of current pathway
    count=0
    countex=0
    countdb=0
    vague=0
    keyReaFound=""
    line=$(echo "$pwyDB" | awk -v i=$i 'NR==i')
    pwy=$(echo "$line" | awk -F "\t" '{print $1}')
    name=$(echo "$line" | awk -F "\t" '{print $2}')
    ecs=$(echo "$line" | awk -F "\t" '{print $7}')
    reaids=$(echo "$line" | awk -F "\t" '{print $6}')
    keyRea=$(echo "$line" | awk -F "\t" '{print $8}' | tr ',' ' ')
    echo -e '\n'$i/$pwyNr: Checking for pathway $pwy $name
    for j in `seq 1 $(echo $ecs | tr "," "\n" | wc -l)`
    #for ec in $(echo $ecs | tr "," "\n")
    do 
        ec=$(echo $ecs | awk -v j=$j -F ',' '{print $j}')
        rea=$(echo $reaids | awk -v j=$j -F ',' '{print $j}')
        re="([0-9]+.[0-9]+.[0-9]+.[0-9]+)"
        test=$(if [[ $ec =~ $re ]]; then echo ${BASH_REMATCH[1]}; fi) # check if not trunked ec number (=> too many hits)
        if [ -n "$test" ]; then
            ((count++))
            getDBhit # get db hits for this reactions
            pwyCandAll="$pwyCandAll$dbhit"
            query=$seqpath$ec.fasta
            if [ ! -s $query ]; then
                python2 $dir/src/uniprot.py "$ec" "$taxonomy" # if sequence data not available then download from uniprot
            fi
            if [ -s $query ]; then
                out=$pwy-$ec.blast
                tblastn -db orgdb -query $query -outfmt '6 qseqid sseqid pident evalue bitscore qcovs stitle' >$out 
                if [ -s $out ]; then
                    bhit=$(cat $out | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff -v covcutoff=$covcutoff '{if ($3>=identcutoff && $5>=bitcutoff && $6>=covcutoff) print $0}')
                    if [ -n "$bhit" ]; then
                        bestIdentity=$(echo "$bhit" | sort -rgk 5,5 | head -1 | cut -f3)
                        bestBitscore=$(echo "$bhit" | sort -rgk 5,5 | head -1 | cut -f5)
                        bestCoverage=$(echo "$bhit" | sort -rgk 5,5 | head -1 | cut -f6)
                        echo -e '\t'Blast hit: $rea $ec "(bit=$bestBitscore, id=$bestIdentity, cov=$bestCoverage)"
                        # check if key reactions of pathway
                        if [[ $keyRea = *"$rea"* ]]; then
                            echo -e '\t\t'Key reaction found!!
                            keyReaFound="$keyReaFound $rea"
                        fi
                        
                        if [ -n "$dbhit" ]; then
                            dbhit="$(echo $dbhit | tr ' ' '\n' | sort | uniq | tr '\n' ' ')" # remove duplicates
                            echo -e '\t\t'Candidate reaction for import: `echo "$dbhit" | wc -w`
                            pwyCand="$pwyCand$dbhit " # remember candidate reaction
                            ((countdb++))
                        else
                            echo -e '\t\t'No candidate reaction found for import
                        fi
                        ((countex++))
                    else
                        someIdentity=$(cat $out | sort -rgk 3,3 | head -1 | cut -f3)
                        someBitscore=$(cat $out | sort -rgk 5,5 | head -1 | cut -f5)
                        someCoverage=$(cat $out | sort -rgk 6,6 | head -1 | cut -f6)
                        echo -e '\t'No significant blast hits found: $rea $ec "\n\t\t(max: id=$someIdentity bit=$someBitscore cov=$someCoverage)"
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
    completness=$(echo "scale=0; 100*$countex/$count" | bc)
    echo "Pathway completness: $countex/$count ($completness%) with $vague reactions of unclear state"
    echo -e Hits with candidate reactions in database: $countdb/$count
    if [ -n "$keyReaFound" ]; then
        CountKeyReaFound=$(echo $keyReaFound | tr ' ' '\n' |  sort | uniq | wc -l)
    else
        CountKeyReaFound=0
    fi
    echo -e Key reactions: $CountKeyReaFound/$(echo $keyRea | wc -w)
    
    if [ $count -ne 0 ] && [ $completness -ge 90 ]; then
        bestCand="$bestCand$pwyCand " # save candidates from almost complete (>=90%) pathways
        bestPwy="$bestPwy$name\n"
    fi
    if [ "$strictCandidates" = false ] && [ $CountKeyReaFound -ge 1 ] && [ $CountKeyReaFound == $(echo $keyRea | wc -w) ] && [ $count -ne 0 ] && [ $completness -ge $completnessCutoff ] && [ $completness -lt 100 ]; then
        echo "Consider pathway to be present because of key enzyme!"
        cand="$cand$pwyCandAll"
        if [[ $bestPwy != *"$name"* ]]; then
           bestCand="$bestCand$pwyCandAll "
           bestPwy="$bestPwy$name ($completness% completness, added because of key enzyme)\n"
        fi
    else
        cand="$cand$pwyCand "
    fi
done

cand="$(echo $cand | tr ' ' '\n' | sort | uniq | tr '\n' ' ')" # remove duplicates
echo -e '\n'Total candidate reactions:
echo $cand

echo -e '\n'Pathways found:
echo -e $bestPwy

bestCand="$(echo $bestCand | tr ' ' '\n' | sort | uniq | tr '\n' ' ')" # remove duplicates
echo -e '\n'Candidate reactions from complete pathways:
echo -e $bestCand

# export found reactions 
#echo $bestCand > newReactions.lst
echo $cand > newReactions.lst
cp newReactions.lst $curdir/${fastaID}-$pathways-Reactions.lst

# add reactions and write new sbml model
if [ -n "$input_sbml" ] ; then # if there is an xml file
    echo ""
    $dir/src/sbml_write.R newReactions.lst $dir
    modelold=$(basename $sbml)
    modelnew="${modelold%.*}G.xml"
    cp modelnew.xml $curdir/$modelnew
fi

