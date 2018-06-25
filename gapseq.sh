#!/bin/bash

# TODO: metacyc superpathways seems to be incomplete e.g. ASPASN-PWY
# TODO: limit pwy search to taxonomic scope
# TODO: save dummy seq file for ec without uniprot hit (save nonsense requests)
# TODO: handle incomplete/unspecific ecs from metacyc (e.g. get ec from kegg, update maually or get genes from metacyc)
# TODO: if taxonomic range is not bacteria, then sequence data must be updated!

pathways=""
database="seed"
verbose=0
taxonomy="Bacteria"
bitcutoff=50 # cutoff blast: min bit score
identcutoff=0   # cutoff blast: min identity
covcutoff=75 # cutoff blast: min coverage
strictCandidates=false
completnessCutoff=66 # consider pathway to be present if other hints (e.g. key enzyme present) are avaiable and pathway completness is at least as high as completnessCutoff (requires strictCandidates=false)
completnessCutoffNoHints=80 # consider pathway to be present if no hints are avaiable (requires stricCandidates=false)
addVague=true # should vague reactions (trunked EC number or no sequence data) be added when there is a hit in reaction DB?
onlyMetacyc=false
blast_format="qseqid pident evalue bitscore qcovs stitle sstart send sseq"
blast_back=false

usage()
{
    echo "Usage"
    echo "$0 -p keyword / -e ec [-d database] [-t taxonomy] file.fasta."
    echo "  -p keywords such as pathways or susbstem (for example amino,nucl,cofactor,carbo,polyamine)"
    echo "  -e search by ec numbers (comma separated)"
    echo "  -d database: vmh or seed (default: $database)"
    echo "  -t taxonomic range (default: $taxonomy)"
    echo "  -b bit score cutoff for local alignment (default: $bitcutoff)"
    echo "  -i identity cutoff for local alignment (default: $identcutoff)"
    echo "  -c coverage cutoff for local alignment (default: $covcutoff)"
    echo "  -s strict candidate reaction handling (do _not_ use pathway completness, key kenzymes and operon structure to infere if imcomplete pathway could be still present (default: $strictCandidates)"
    echo "  -n Add vague reactions (i.e. EC is trunked and no sequences is available) if pathway is otherwise complete (default: $addVague)"
    echo "  -o use only MetaCyc pathway database (default: $onlyMetacyc)"
    echo "  -u suffix used for output files (default: pathway keyword)"
    echo "  -a blast hits back against uniprot enzyme database"
exit 1
}


# paths and variables
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
uniprotIdentity=0.9 # clustered uniprot database (0.5 or 0.9)
metaPwy=$dir/dat/meta_pwy.tbl
keggPwy=$dir/dat/kegg_pwy.tbl
metaRea=$dir/dat/meta_rea.tbl
reaDB1=$dir/dat/vmh_reactions.csv
reaDB2=$dir/dat/bigg_reactions.tbl
reaDB3=$dir/dat/seed_reactions.tsv
reaDB4=$dir/dat/mnxref_seed.tsv
brenda=$dir/dat/brenda_ec.csv
seedEC=$dir/dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv


# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "h?p:e:d:i:b:c:vs:t:snou:a" opt; do
    case "$opt" in
    h|\?)
        usage
        exit 0
        ;;
    p)  
        pathways=$OPTARG
        ;;
    e)  
        ecnumber=$OPTARG
        ;;
    d)  
        database=$OPTARG
        ;;
    v)  
        verbose=1
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
    s)
        strictCandidates=true
        ;;
    n)
        addVague=false
        ;;
    o)
        onlyMetacyc=true
        ;;
    u)
        output_suffix=$OPTARG
        ;;
    a)
        blast_back=true
        ;;
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

# pathways or ec number as well as fasta file have to be provided
( [ -z "$pathways" ] && [ -z "$ecnumber" ] ) || [ -z "$fasta" ]  && { usage; }

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
        pwyKey="Amino-Acid-Biosynthesis|Nucleotide-Biosynthesis|Cofactor-Biosynthesis|Carbohydrates-Degradation|CARBO-BIOSYNTHESIS|Polyamine-Biosynthesis|Fatty-acid-biosynthesis|Energy-Metabolism|Terpenoid-Biosynthesis|Chorismate-Biosynthesis"
        ;;
    kegg)
        pwyKey=kegg
        ;;
    *)
        pwyKey=$pathways
        ;;
esac


# squence directory
seqpath=$dir/dat/seq/$taxonomy/unipac$(printf %.0f $(echo "$uniprotIdentity * 100" | bc -l))
mkdir -p $seqpath

# tmp working directory
cd $(mktemp -d)


if [ -n "$ecnumber" ]; then
    # create dummpy pwy template for given ec number
    pwyKey=$ecnumber
    pwyDB=$(echo -e "dummy\t$ecnumber\t\t\t\t$ecnumber\t$ecnumber")
    pathways="ec"
else
    # get entries for pathways from databases
    if [ "$onlyMetacyc" == true ]; then
        cat $metaPwy > allPwy
    else
        cat $metaPwy $keggPwy > allPwy
    fi
    pwyDB=$(cat allPwy | grep -wEi $pwyKey)
    [ -z "$ecnumber" ] && [ -z "$pwyDB" ] && { echo "No pathways found for key $pwyKey"; exit 1; }
fi
[ -z "$output_suffix" ] && output_suffix=$pathways


# function to get database hits for ec number
getDBhit(){
    kegg=$(grep -wF $rea $metaRea | awk -F "\t" {'print $5'})
    altec=$(grep -wF $ec $brenda | grep -P "([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)" -o | grep -v $ec)
    
    # 1) search in reaction db by EC
    if [ "$database" == "vmh" ]; then
        dbhit=$(grep -wF $ec $reaDB1 | awk -F ',' '{print $1}')
    elif [ "$database" == "seed" ]; then
        dbhit=$(grep -wF $ec $seedEC | awk -F '\t' '{print $1}' | tr '|' ' ')
        dbhit="$dbhit $(grep -wF $ec $reaDB4 | awk -F '\t' '{print $4}' | tr '\n' ' ')"
    fi

    # 2) search in reaction db by kegg identifier 
    if [ "$database" == "vmh" ]; then
        [ -n "$kegg" ]  && dbhit="$dbhit $(grep -wE "$(echo $kegg |tr ' ' '|')" $reaDB1 | awk -F ',' '{print $1}')"
    elif [ "$database" == "seed" ]; then
        [ -n "$kegg" ] && dbhit="$dbhit $(grep -wE "$(echo $kegg | tr ' ' '|')" $reaDB3 | awk -F '\t' '$18 == "OK" {print $1}' )" # only consider reactions which are OK
        [ -n "$kegg" ] && dbhit="$dbhit $(grep -wE "$(echo $kegg | tr ' ' '|')" $reaDB4 | awk -F '\t' '{print $4}' | tr '\n' ' ')"
    fi

    # 3) search in reaction db by alternative EC
    if [ "$database" == "vmh" ]; then
        [ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo $altec | tr ' ' '|')" $reaDB1 | awk -F ',' '{print $1}')" # take care of multiple EC numbers
    elif [ "$database" == "seed" ]; then
        [ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo $altec | tr ' ' '|')" $seedEC | awk -F '\t' '{print $1}' | tr '|' ' ')" # take care of multiple EC numbers
        [ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo $altec | tr ' ' '|')" $reaDB4 | awk -F '\t' '{print $4}')" # take care of multiple EC numbers
    fi
    
    # 4) search in bigg db by metacyc id (does only make sense for vmh/bigg namespace)
    if [ "$database" == "vmh" ]; then
        dbhit="$dbhit $(grep -wF $rea $reaDB2 | awk '{print $1}')"
    fi

    [ "$dbhit" == " " ] && dbhit=""
}





# create blast database
makeblastdb -in $fasta -dbtype nucl -out orgdb >/dev/null


cand=""     #list of candidate reactions to be added
bestCand="" # list of candidates from (almost) complete pathways
bestPwy=""  # list of found pathways
echo -e "ID\tName\tPrediction\tCompletness\tVagueReactions\tKeyReactions\tKeyReactionsFound\tReactionsFound" > output.tbl # pahtway statistics file

pwyNr=$(echo "$pwyDB" | wc -l)
echo Checking for pathways and reactions in: $1 $pwyKey
echo Number of pathways to be considered: $pwyNr
for i in `seq 1 $pwyNr`
do
    pwyCand="" # candidate reaction of current pathway
    pwyCandAll="" # all possible reaction of current pathway
    pwyVage="" # reaction belonging to trunked EC numbers (no sequence blast possible..)
    count=0
    countex=0
    countexList="" # list with reactions ids found
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
        re="([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)"
        test=$(if [[ $ec =~ $re ]]; then echo ${BASH_REMATCH[1]}; fi) # check if not trunked ec number (=> too many hits)
        ((count++))
        if [ -n "$test" ]; then
            getDBhit # get db hits for this reactions
            pwyCandAll="$pwyCandAll$dbhit "
            query=$seqpath/$ec.fasta
            if [ ! -f $query ]; then # check if sequence is not available => try to download
                #python2 $dir/src/uniprot.py "$ec" "$taxonomy" # if sequence data not available then download from uniprot
                
                echo -e '\t'Downloading sequence for: $ec 
                $dir/src/uniprot.sh -e "$ec" -t "$taxonomy" -i $uniprotIdentity >/dev/null
            fi
            if [ -s $query ]; then
                out=$ec.blast
                if [ ! -f $out ]; then # check if there is a former hit
                    csplit -s -z $query '/>/' '{*}' # split multiple fasta file and skip further testing if a significant hit is found
                    for q in `ls xx*`
                    do
                        tblastn -db orgdb -query $q -outfmt "6 $blast_format" >>$out 
                        bhit=$(cat $out | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff -v covcutoff=$covcutoff '{if ($2>=identcutoff && $4>=bitcutoff && $5>=covcutoff) print $0}')
                        if [ -n "$bhit" ]; then
                            break
                        fi
                    done
                    rm xx*
                fi
                if [ -s $out ]; then
                    bhit=$(cat $out | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff -v covcutoff=$covcutoff '{if ($2>=identcutoff && $4>=bitcutoff && $5>=covcutoff) print $0}')
                    if [ -n "$bhit" ]; then
                        bestIdentity=$(echo "$bhit" | sort -rgk 4,4 | head -1 | cut -f2)
                        bestBitscore=$(echo "$bhit" | sort -rgk 4,4 | head -1 | cut -f4)
                        bestCoverage=$(echo "$bhit" | sort -rgk 4,4 | head -1 | cut -f5)
                        besthit_all=$(echo "$bhit" | sort -rgk 4,4 | head -1)
                        echo -e "$rea\t$ec\t$besthit_all" >> reactions.tbl 
                        echo -e '\t'Blast hit: $rea $ec "(bit=$bestBitscore, id=$bestIdentity, cov=$bestCoverage)"
                        # check if key reactions of pathway
                        if [[ $keyRea = *"$rea"* ]]; then
                            echo -e '\t\t'Key reaction found!!
                            keyReaFound="$keyReaFound $rea"
                        fi
                        #blast hit back to uniprot enzyme database
                        if [ "$blast_back" = true ]; then
                            echo "$bhit" | sort -rgk 4,4 | head -1 | cut -f9 | sed 's/-/*/g' > "$ec.hit.fasta"
                            echo "Blast best hit against uniprot db:"
                            blastp -db $dir/dat/seq/uniprot_sprot -query "$ec.hit.fasta" -outfmt '6 pident bitscore qcovs stitle' | awk '{if ( $4>50 ) print $0}' | sort -rgk 2,2 | head -n 3
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
                        countexList="$countexList$rea "
                    else
                        someIdentity=$(cat $out | sort -rgk 4,4 | head -1 | cut -f2)
                        someBitscore=$(cat $out | sort -rgk 4,4 | head -1 | cut -f4)
                        someCoverage=$(cat $out | sort -rgk 4,4 | head -1 | cut -f5)
                        somehit_all=$( cat $out | sort -rgk 4,4 | head -1)
                        echo -e "$rea\t$ec\t$somehit_all" >> reactions.tbl 
                        echo -e '\t'No significant blast hits found: $rea $ec "\n\t\t(best one: id=$someIdentity bit=$someBitscore cov=$someCoverage)"
                        if [[ -n "$dbhit" ]];then
                            dbhit="$(echo $dbhit | tr ' ' '\n' | sort | uniq | tr '\n' ' ')" # remove duplicates
                            pwyNoSeqFound="$pwyNoSeqFound$dbhit "
                        fi
                    fi
                else
                    echo -e '\t'No blast hits found for reaction: $rea $ec
                    getDBhit
                    if [[ -n "$dbhit" ]];then
                        dbhit="$(echo $dbhit | tr ' ' '\n' | sort | uniq | tr '\n' ' ')" # remove duplicates
                        pwyNoSeqFound="$pwyNoSeqFound$dbhit "
                    fi
                fi
            else
                echo -e '\t'No sequence data found for $rea $ec ..skipping..
                ((vague++))
                getDBhit
                if [[ -n "$dbhit" ]];then
                    dbhit="$(echo $dbhit | tr ' ' '\n' | sort | uniq | tr '\n' ' ')" # remove duplicates
                    pwyNoSeqFound="$pwyNoSeqFound$dbhit "
                fi
            fi
        else
            echo -e '\t'EC number too unspecific: $rea $ec ..skipping..
            ((vague++))
            if [ -n "$ec" ]; then
                ec=$ec.-
                getDBhit
                [[ -n "$dbhit" ]] && pwyVage="$pwyVage$dbhit "
            fi
        fi
        [[ verbose -gt 0 ]] && echo -e "\t\tCandidate reactions: $dbhit"
    done # pathway
    
    if [ $count -eq 0 ]; then
        completness=0
    else
        check_vague=$(echo "$vague < $count*0.5" | bc) # vague reactions shouldn't make more than half of total reactions
        if [ "$addVague" = true ] && [ $check_vague -eq 1 ] ; then # if vague reaction are considered they should not influence the completness treshold
            completness=$(echo "scale=0; 100*($countex+$vague)/$count" | bc)
        else
            completness=$(echo "scale=0; 100*($countex)/$count" | bc)
        fi
    fi
    if [ $vague -eq 0 ] || [ "$addVague" = false ] || [ $check_vague -eq 0 ]; then
        echo "Pathway completness: $countex/$count ($completness%)"
    else
        echo "Pathway completness: ($countex+$vague)/$count ($completness%) with $vague reactions of unclear state"
    fi

    echo -e Hits with candidate reactions in database: $countdb/$count
    if [ -n "$keyReaFound" ]; then
        CountKeyReaFound=$(echo $keyReaFound | tr ' ' '\n' |  sort | uniq | wc -l)
    else
        CountKeyReaFound=0
    fi
    CountKeyRea=$(echo $keyRea | wc -w)
    CountTotalKeyRea=$(echo $keyRea | wc -w)
    echo -e Key reactions: $CountKeyReaFound/$CountKeyRea
    
    # add reactions of pathways (even if no blast hit) if above trashold (and no key enzyme is missed)
    prediction=false
    if [ $count -ne 0 ] && [ $completness -ge $completnessCutoffNoHints ] && [ $CountKeyReaFound -eq $CountTotalKeyRea ]; then
        echo "Consider pathway to be present because of completness treshold!"
        prediction=true
        bestCand="$bestCand$pwyCand " # save candidates from almost complete pathways
        if [[ -n "$pwyVage" ]] && [[ "$addVague" = true ]]; then
           bestCand="$bestCand$pwyVage$pwyNoSeqFound " # add vague reaction for pathways that are present
           cand="$cand$pwyVage$pwyNoSeqFound "
        fi
        if [[ -n "$pwyNoSeqFound" ]] && [[ "$strictCandidates" = false ]]; then
           bestCand="$bestCand$pwyNoSeqFound " # add reaction with no blast hits
           cand="$cand$pwyNoSeqFound "
        fi
        if [[ $completness -lt 100 ]]; then
            bestPwy="$bestPwy$name ($completness% completness, added because of treshhold)\n"
        else
            bestPwy="$bestPwy$name\n"
        fi
    fi
    if [ $CountKeyReaFound -ge 1 ] && [ $CountKeyReaFound -eq $CountTotalKeyRea ] && [ $count -ne 0 ] && [ $completness -ge $completnessCutoff ] && [ $completness -lt 100 ]; then
        echo "Consider pathway to be present because of key enzyme!"
        prediction=true
        cand="$cand$pwyCandAll"
        if [[ $bestPwy != *"$name"* ]]; then # if not alrady added because of completness (s.a.)
           bestCand="$bestCand$pwyCandAll "
           if [[ -n "$pwyVage" ]] && [[ "$addVague" = true ]]; then
              bestCand="$bestCand$pwyVage " # add vague reaction for pathways that are present
              cand="$cand$pwyVage "
            fi
            if [[ -n "$pwyNoSeqFound" ]] && [[ "$strictCandidates" = false ]]; then
               bestCand="$bestCand$pwyNoSeqFound " # add reaction with no blast hits
               cand="$cand$pwyNoSeqFound "
            fi
           bestPwy="$bestPwy$name ($completness% completness, added because of key enzyme)\n"
        fi
    else
        cand="$cand$pwyCand "
    fi
    
    echo -e "$pwy\t$name\t$prediction\t$completness\t$vague\t$CountKeyRea\t$CountKeyReaFound\t$countexList" >> output.tbl # write down some statistics

done

cand="$(echo $cand | tr ' ' '\n' | sort | uniq | tr '\n' ' ')" # remove duplicates
if [[ verbose -gt 0 ]]; then
    echo -e '\n'Total candidate reactions:
    echo $cand
fi

echo -e '\n'Pathways found:
echo -e $bestPwy

bestCand="$(echo $bestCand | tr ' ' '\n' | sort | uniq | tr '\n' ' ')" # remove duplicates
if [[ verbose -gt 0 ]]; then
    echo -e '\n'Candidate reactions from complete pathways:
    echo -e $bestCand
fi

# export found reactions 
echo -e Candidate reactions found: $(echo "$cand" | wc -w) '\n'
#echo $bestCand > newReactions.lst
echo $cand > newReactions.lst
cp newReactions.lst $curdir/${fastaID}-$output_suffix-Reactions.lst
cp output.tbl $curdir/${fastaID}-$output_suffix-Pathways.tbl
[ -f reactions.tbl ] && echo "#rxn ec $blast_format" | cat - reactions.tbl | awk '!a[$0]++' > $curdir/${fastaID}-$output_suffix-blast.tbl # add header and remove duplicates


ps -p $$ -o %cpu,%mem,cmd
times
