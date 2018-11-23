#!/bin/bash

# TODO: metacyc superpathways seems to be incomplete e.g. ASPASN-PWY
# TODO: limit pwy search to taxonomic scope
# TODO: save dummy seq file for ec without uniprot hit (save nonsense requests)
# TODO: handle incomplete/unspecific ecs from metacyc (e.g. get ec from kegg, update maually or get genes from metacyc)
# TODO: if taxonomic range is not bacteria, then sequence data must be updated!

start_time=`date +%s`

pathways=""
database="seed"
pwyDatabase="metacyc,custom"
verbose=1
taxonomy="Bacteria"
bitcutoff=50 # cutoff blast: min bit score
identcutoff=0   # cutoff blast: min identity
identcutoff_exception=70  # min identity for enzymes marked as false friends (hight seq similarity but different function)
covcutoff=75 # cutoff blast: min coverage
subunit_cutoff=50 # more than this % of subunits must be found 
strictCandidates=false
completenessCutoff=66 # consider pathway to be present if other hints (e.g. key enzyme present) are avaiable and pathway completeness is at least as high as completenessCutoff (requires strictCandidates=false)
completenessCutoffNoHints=80 # consider pathway to be present if no hints are avaiable (requires stricCandidates=false)
blast_back=false
noSuperpathways=true
vagueCutoff=0.3 # cutoff for vague reactions. If the amount of vague reactions in a pathways is more then this their influence will not be recognized even with strictCandidates=false
onlyList=false
skipBlast=false
includeSeq=false

usage()
{
    echo "Usage"
    echo "$0 -p keyword / -e ec [-d database] [-t taxonomy] file.fasta."
    echo "  -p keywords such as pathways or subsystems (for example amino,nucl,cofactor,carbo,polyamine)"
    echo "  -e search by ec numbers (comma separated)"
    echo "  -r search by enzyme name (colon separated)"
    echo "  -d database: vmh or seed (default: $database)"
    echo "  -t taxonomic range for sequences to be downloaded (default: $taxonomy)"
    echo "  -b bit score cutoff for local alignment (default: $bitcutoff)"
    echo "  -i identity cutoff for local alignment (default: $identcutoff)"
    echo "  -c coverage cutoff for local alignment (default: $covcutoff)"
    echo "  -s strict candidate reaction handling (do _not_ use pathway completeness, key kenzymes and operon structure to infere if imcomplete pathway could be still present (default: $strictCandidates)"
    echo "  -u suffix used for output files (default: pathway keyword)"
    echo "  -a blast hits back against uniprot enzyme database"
    echo "  -n Consider superpathways of metacyc database"
    echo "  -l Select the pathway database (MetaCyc, KEGG, SEED, all; default: $pwyDatabase)"
    echo "  -o Only list pathways found for keyword; default $onlyList)"
    echo "  -x Do not blast only list pathways, reactions and check for available sequences; default $skipBlast"
    echo "  -q Include sequences of hits in log files; default $includeSeq"

    echo "  -v verbose level, 0 for nothing, 1 for pathway infos, 2 for full (default $verbose)"
exit 1
}


# paths and variables
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
uniprotIdentity=0.9 # clustered uniprot database (0.5 or 0.9)
metaPwy=$dir/dat/meta_pwy.tbl
keggPwy=$dir/dat/kegg_pwy.tbl
seedPwy=$dir/dat/seed_pwy.tbl
customPwy=$dir/dat/custom_pwy.tbl
metaRea=$dir/dat/meta_rea.tbl
reaDB1=$dir/dat/vmh_reactions.csv
reaDB2=$dir/dat/bigg_reactions.tbl
reaDB3=$dir/dat/seed_reactions.tsv
reaDB4=$dir/dat/mnxref_seed.tsv
reaDB5=$dir/dat/mnxref_seed-other.tsv
brenda=$dir/dat/brenda_ec.csv
seedEC=$dir/dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv
seedEnzymesNames=$dir/dat/seed_Enzyme_Name_Reactions_Aliases.tsv


# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "h?p:e:r:d:i:b:c:v:st:snou:al:oxq" opt; do
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
    r)  
        reaname=$OPTARG
        ;;
    d)  
        database=$OPTARG
        ;;
    v)  
        verbose=$OPTARG
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
    u)
        output_suffix=$OPTARG
        ;;
    a)
        blast_back=true
        ;;
    n)
        noSuperpathways=false
        ;;
    l)
        pwyDatabase=$OPTARG
        ;;
    o)
        onlyList=true
        ;;
    x)
        skipBlast=true
        ;;
    q)
        includeSeq=true
        ;;
    esac
done
shift $((OPTIND-1))
[ "$1" = "--" ] && shift

# after parsing arguments, only fasta file shoud be there
[ "$#" -ne 1 ] && { usage; }

# blast format
if [ "$includeSeq" = true ]; then
    blast_format="qseqid pident evalue bitscore qcovs stitle sstart send sseq"
else
    blast_format="qseqid pident evalue bitscore qcovs stitle sstart send"
fi

# tmp working directory
fasta=$(readlink -f "$1") # save input file before changing to temporary directory
cd $(mktemp -d)

# get fasta file
if [[ $fasta == *.gz ]]; then # in case fasta is in a archive
    tmp_fasta=$(basename "${fasta}" .gz)
    gunzip -c $fasta > $tmp_fasta
    fasta=$tmp_fasta
fi
[[ ! -s $fasta ]] && { echo Invalid file: $1; exit 0; }
tmpvar=$(basename $fasta)
fastaID="${tmpvar%.*}"

# pathways or ec number as well as fasta file have to be provided
( [ -z "$pathways" ] && [ -z "$ecnumber" ] && [ -z "$reaname" ] ) || [ -z "$fasta" ]  && { usage; }

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


if [ -n "$ecnumber" ] || [ -n "$reaname" ]; then
    # create dummpy pwy template for given ec number
    if [[ -z "$ecnumber" ]]; then
        rea_count=$(echo $reaname | tr ';' '\n' | wc -l)
        ecnumber=$(echo $reaname | grep -o ";" | tr ';' ',') # get dummy empty comma seperated ec numbers
    elif [[ -z "$reaname" ]]; then
        rea_count=$(echo $ecnumber | tr ',' '\n' | wc -l)
        reaname=$(echo $ecnumber | grep -o "," | tr ',' ';') # get dummy empty colon seperated reaction names
    else  
        rea_count=$(echo $ecnumber | tr ',' '\n' | wc -l)
    fi
    rea_id=$(seq 1 $rea_count | awk '{print "reaction"$1}' |tr '\n' ',')
    pwyKey="custom"
    pwyDB=$(echo -e "custom\t$ecnumber\t\t\t\t${rea_id::-1}\t$ecnumber\t\t$reaname")
    pathways="custom"
else
    pwyDatabase=$(echo $pwyDatabase | tr '[:upper:]' '[:lower:]')
    # get entries for pathways from databases
    [[ verbose -ge 1 ]] && { echo $pwyDatabase; }
    [[ "$pwyDatabase" =~ "all" ]]     && cat $metaPwy $keggPwy $seedPwy $customPwy > allPwy
    [[ "$pwyDatabase" =~ "metacyc" ]] && cat $metaPwy >> allPwy
    [[ "$pwyDatabase" =~ "kegg" ]]    && cat $keggPwy >> allPwy
    [[ "$pwyDatabase" =~ "seed" ]]    && cat $seedPwy >> allPwy
    [[ "$pwyDatabase" =~ "custom" ]]  && cat $customPwy >> allPwy
    pwyDB=$(cat allPwy | grep -wEi $pwyKey)
    [[ "$noSuperpathways" = true ]] && pwyDB=$(echo "$pwyDB" | grep -v 'Super-Pathways')
    [ -z "$ecnumber" ] && [ -z "$pwyDB" ] && { echo "No pathways found for key $pwyKey"; exit 1; }
fi
[ -z "$output_suffix" ] && output_suffix=$pathways


# function to get database hits for ec number
getDBhit(){
    kegg=$(grep -wFe "$rea" $metaRea | awk -F "\t" {'print $5'})

    # 1) search in reaction db by EC
    if [[ -n "$EC_test" ]]; then
        if [ "$database" == "vmh" ]; then
            dbhit=$(grep -wF $ec $reaDB1 | awk -F ',' '{print $1}')
        elif [ "$database" == "seed" ]; then
            dbhit=$(grep -wF $ec $seedEC | awk -F '\t' '{print $1}' | tr '|' ' ')
            dbhit="$dbhit $(grep -wF $ec $reaDB4 | awk -F '\t' '{print $4}' | tr '\n' ' ')"
        fi
    fi

    # 2) search in reaction db by kegg identifier 
    if [ "$database" == "vmh" ]; then
        [ -n "$kegg" ]  && dbhit="$dbhit $(grep -wE "$(echo $kegg |tr ' ' '|')" $reaDB1 | awk -F ',' '{print $1}')"
    elif [ "$database" == "seed" ]; then
        [ -n "$kegg" ] && dbhit="$dbhit $(grep -wE "$(echo $kegg | tr ' ' '|')" $reaDB3 | awk -F '\t' '$18 == "OK" {print $1}' )" # only consider reactions which are OK
        [ -n "$kegg" ] && dbhit="$dbhit $(grep -wE "$(echo $kegg | tr ' ' '|')" $reaDB4 | awk -F '\t' '{print $4}' | tr '\n' ' ')"
    fi

    # 3) search in reaction db by alternative EC
    if [[ -n "$EC_test" ]]; then
        altec=$(grep -wF $ec $brenda | grep -P "([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)" -o | grep -v $ec)
        if [ "$database" == "vmh" ]; then
            [ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo $altec | tr ' ' '|')" $reaDB1 | awk -F ',' '{print $1}')" # take care of multiple EC numbers
        elif [ "$database" == "seed" ]; then
            [ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo $altec | tr ' ' '|')" $seedEC | awk -F '\t' '{print $1}' | tr '|' ' ')" # take care of multiple EC numbers
            [ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo $altec | tr ' ' '|')" $reaDB4 | awk -F '\t' '{print $4}')" # take care of multiple EC numbers
        fi
    fi

    # 4) search in bigg db by metacyc id (does only make sense for vmh/bigg namespace)
    if [ "$database" == "vmh" ]; then
        dbhit="$dbhit $(grep -wFe "$rea" $reaDB2 | awk '{print $1}')"
    fi

    # 5) match reaction using mnxref namespace
    if [ "$database" == "seed" ]; then
        dbhit="$dbhit $(grep -wFe "$rea" $reaDB5 | awk '{print $2}')"
    fi
   
    # 6) match reaction using custom enzyme-name - seedID mapping
    if [ "$database" == "seed" ]; then
        dbhit="$dbhit $(grep -wFe "$reaName" $seedEnzymesNames | awk -F '\t' ' {print $1}')"
    fi

    [ "$dbhit" == " " ] && dbhit=""
}





# create blast database
makeblastdb -in $fasta -dbtype nucl -out orgdb >/dev/null


cand=""     #list of candidate reactions to be added
bestPwy=""  # list of found pathways
echo -e "ID\tName\tPrediction\tCompleteness\tVagueReactions\tKeyReactions\tKeyReactionsFound\tReactionsFound" > output.tbl # pahtway statistics file

pwyNr=$(echo "$pwyDB" | wc -l)
[[ verbose -ge 1 ]] && echo Checking for pathways and reactions in: $1 $pwyKey
[[ verbose -ge 1 ]] && echo Number of pathways to be considered: $pwyNr
for i in `seq 1 $pwyNr`
do
    pwyCand="" # candidate reaction of current pathway
    pwyVage="" # reaction belonging to trunked EC numbers (no sequence blast possible..)
    pwyNoHitFound="" # remember reactions without blast hit so that they can be added in case of high pathway completeness 
    count=0
    countex=0
    countexList="" # list with reactions ids found
    countdb=0
    vague=0
    keyReaFound=""
    vagueKeyReaFound=""
    line=$(echo "$pwyDB" | awk -v i=$i 'NR==i')
    pwy=$(echo "$line" | awk -F "\t" '{print $1}')
    name=$(echo "$line" | awk -F "\t" '{print $2}')
    ecs=$(echo "$line" | awk -F "\t" '{print $7}')
    reaids=$(echo "$line" | awk -F "\t" '{print $6}')
    reaNames=$(echo -e "$line" | awk -F "\t" '{print $9}')
    keyRea=$(echo "$line" | awk -F "\t" '{print $8}' | tr ',' ' ')
    pwyHierarchy=$(echo "$line" | awk -F "\t" '{print $4}' | sed 's/|\|THINGS\|Generalized-Reactions\|Pathways\|FRAMES//g' | sed 's/,,//g')
    [[ verbose -ge 1 ]] && echo -e '\n'$i/$pwyNr: Checking for pathway $pwy $name
    [[ verbose -ge 1 ]] && echo "($pwyHierarchy)"
    [[ "$onlyList" = true ]] && { continue; }
    for j in `seq 1 $(echo $ecs | tr "," "\n" | wc -l)`
    #for ec in $(echo $ecs | tr "," "\n")
    do 
        dbhit=""
        ec=$(echo $ecs | awk -v j=$j -F ',' '{print $j}')
        rea=$(echo $reaids | awk -v j=$j -F ',' '{print $j}')
        reaName=$(echo $reaNames | awk -v j=$j -F ';' '{print $j}' | tr -d '|')
        re="([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)"
        EC_test=$(if [[ $ec =~ $re ]]; then echo ${BASH_REMATCH[1]}; fi) # check if not trunked ec number (=> too many hits)
        [[ verbose -ge 1 ]] && echo -e '\t'$rea $reaName $ec
        [[ -z "$rea" ]] && { continue; }
        [[ -n "$ec" ]] && [[ -n "$reaName" ]] && [[ -n "$EC_test" ]] && { is_exception=$(grep -Fw -e "$ec" -e "$reaName" $dir/dat/exception.tbl | wc -l); }
        ( [[ -z "$ec" ]] || [[ -z "$EC_test" ]] ) && [[ -n "$reaName" ]] && { is_exception=$(grep -Fw "$reaName" $dir/dat/exception.tbl | wc -l); }
        [[ -n "$ec" ]] && [[ -z "$reaName" ]] && [[ -n "$EC_test" ]] && { is_exception=$(grep -Fw "$ec" $dir/dat/exception.tbl | wc -l); }
        if [[ $is_exception -gt 0 ]] && [[ $identcutoff -lt $identcutoff_exception ]];then # take care of similair enzymes with different function
            identcutoff_tmp=$identcutoff_exception
            [[ verbose -ge 1 ]] && echo -e "\t\tUsing higher identity cutoff for $rea"
        else
            identcutoff_tmp=$identcutoff
        fi
        ((count++))
        getDBhit # get db hits for this reactions
        dbhit="$(echo $dbhit | tr ' ' '\n' | sort | uniq | tr '\n' ' ')" # remove duplicates
        if [[ -n "$EC_test" ]]; then
            query=$seqpath/$ec.fasta
            if [ ! -f "$query" ]; then # check if sequence is not available => try to download
                [[ verbose -ge 1 ]] && echo -e '\t'Downloading sequence for: $ec 
                $dir/src/uniprot.sh -e "$ec" -t "$taxonomy" -i $uniprotIdentity >/dev/null
            fi
        fi
        # if no EC number is available or no sequence was found for EC number then use reaction name instead for sequence search
        if [[ -n "$reaName" ]] && ( [[ -z "$EC_test" ]] || [[ ! -s "$query" ]] );then
            reaNameHash=$(echo -n "$reaName" | md5sum | awk '{print $1}')
            query="$seqpath/$reaNameHash.fasta"
            if [ ! -f "$query" ]; then # check if sequence is not available => try to download
                [[ verbose -ge 1 ]] && echo -e '\t'Downloading sequence for: $reaName "\n\t\t(hash: $reaNameHash)" 
                $dir/src/uniprot.sh -r "$reaName" -t "$taxonomy" -i $uniprotIdentity >/dev/null
            fi
        fi
        [[ "$skipBlast" = true ]] && { echo -e "\t"$rea $reaName $ec; continue; }

        if [ -s "$query" ]; then
            [[ verbose -ge 1 ]] && echo -e "\t\t$query" 
            #out=$rea.blast #$ec.blast
            query_id=$(basename $query)
            out="${query_id%.*}".blast
            subunits_found=0
            if [ ! -f $out ]; then # check if there is a former hit
                #subunits=$(cat $query | sed -n 's/^>//p' | grep -oP 'subunit [0-9]' | sort | uniq) # check for subunits
                subunits=$(cat $query | sed -n 's/^>//p' | grep -oE 'subunit [0-9]|(alpha|beta|gamma|delta|epsilon) subunit' | sort | uniq) # check for subunits
                subunits=$(echo -e "$subunits\nother" | sed '/^$/d') # add default
                iterations=$(echo -e "$subunits"| wc -l) # every subunit will get a own iteration
                [[ $iteractions -gt 1 ]] && echo -e '\t\t'check subunits: $subunits
                for iter in `seq $iterations`
                do
                    if [ $iterations -gt 1 ]; then 
                        # apt install exonerate
                        fastaindex $query query.idx
                        if [ "$subunit_key" == "other" ];then
                            subunit_key=$(echo "$subunits" | tr '\n' '|' | sed 's/|$//g')
                            cat $query | sed -n 's/^>//p' | grep -vE "$subunit_key" | awk '{print $1}' | sed 's/^>//g' > query_subunit_header
                        else
                            subunit_key=$(echo "$subunits" | sed -n ${iter}p)
                            cat $query | grep "$subunit_key" | awk '{print $1}' | sed 's/^>//g' > query_subunit_header
                        fi
                        #echo -e $iter "\n"
                        fastafetch -f $query -i query.idx -Fq <(sort -u query_subunit_header) > query_subunit.fasta
                        rm query.idx query_subunit_header
                    else
                        cp $query query_subunit.fasta
                    fi
                    #csplit -s -z query_subunit.fasta '/>/' '{*}' # split multiple fasta file and skip further testing if a significant hit is found
                    $dir/src/fasta-splitter.pl --n-parts 10 query_subunit.fasta >/dev/null
                    #for q in `ls xx*`
                    for q in `ls query_subunit.part-*.fasta`
                    do
                        if ! [ -x "$(command -v parallel)" ]; then # try to use parallelized version
                            tblastn -db orgdb -query $q -qcov_hsp_perc $covcutoff -outfmt "6 $blast_format" > query.blast
                        else
                            cat $q | parallel --will-cite --block 50k --recstart '>' --pipe tblastn -db orgdb -qcov_hsp_perc $covcutoff -outfmt \'"6 $blast_format"\' -query - > query.blast
                        fi
                        cat query.blast >> $out
                        bhit=$(cat query.blast | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff_tmp -v covcutoff=$covcutoff '{if ($2>=identcutoff && $4>=bitcutoff && $5>=covcutoff) print $0}')
                        rm query.blast
                        if [ -n "$bhit" ]; then
                            bestsubunithit=$(echo "$bhit" | sort -rgk 4,4 | head -1)
                            hit_id=$(echo $bestsubunithit | awk '{print $1}')
                            #echo -e '\t'subunit $iter found: $hit_id
                            [[ $iterations -gt 1 ]] && cat $q | head -1 | sed "s/^/\t\t$subunit_key hit: /" 
                            ((subunits_found++))
                            break
                        fi
                    done
                    rm query_subunit.part-*.fasta*
                done
                [[ $iterations -gt 1 ]] && [[ verbose -ge 1 ]] &&  echo -e '\t\t'total subunits found: $subunits_found / $iterations
                echo -e $out'\t'$subunits_found'\t'$iterations >> subunits.log # save subunits found
            else
                # get subunit fraction from former run
                subunits_found=$(cat subunits.log | awk -F "\t" -v out=$out '$1==out { print $2 }')
                iterations=$(cat subunits.log | awk -F "\t" -v out=$out '$1==out { print $3 }')
                #echo test:$subunits_found $iterations
            fi
            if [ -s $out ]; then
                bhit=$(cat $out | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff_tmp -v covcutoff=$covcutoff '{if ($2>=identcutoff && $4>=bitcutoff && $5>=covcutoff) print $0}')
                subunit_fraction=$(echo "100*$subunits_found/$iterations" | bc)
                if [ -n "$bhit" ] && [ $subunit_fraction -gt $subunit_cutoff ] ; then
                    bestIdentity=$(echo "$bhit" | sort -rgk 4,4 | head -1 | cut -f2)
                    bestBitscore=$(echo "$bhit" | sort -rgk 4,4 | head -1 | cut -f4)
                    bestCoverage=$(echo "$bhit" | sort -rgk 4,4 | head -1 | cut -f5)
                    besthit_all=$(echo "$bhit" | sort -rgk 4,4 | head -3)
                    bhit_count=$(echo "$bhit" | wc -l)
                    [[ verbose -ge 1 ]] && echo -e '\t\t'Blast hit \(${bhit_count}x\)
                    [[ verbose -ge 1 ]] && echo "$besthit_all" | awk '{print "\t\t\tbit="$4 " id="$2 " cov="$5 " hit="$1}'
                    # check if key reactions of pathway
                    is_bidihit=NA
                    if [[ $keyRea = *"$rea"* ]]; then
                        [[ verbose -ge 1 ]] && echo -e '\t\t--> KEY reaction found <--'
                        keyReaFound="$keyReaFound $rea"
                    fi
                    #blast hit back to uniprot enzyme database
                    if [ "$blast_back" = true ]; then
                        echo "$bhit" | sort -rgk 4,4 | head -3 | cut -f9 | sed 's/-/*/g' > "$rea.hit.fasta"
                        [[ verbose -ge 1 ]] && echo -e "\t\tBlast best hit against uniprot db:"
                        blastp -db $dir/dat/seq/uniprot_sprot -query "$rea.hit.fasta" -outfmt '6 pident bitscore qcovs sseqid qseqid' > $rea.hit.blast 

                        forward_hit=$(echo "$bhit" | sort -rgk 4,4 | cut -f1 | sed 's/UniRef90_//g' | sort | uniq | tr '\n' '|' | sed 's/|$//g')

                        back_hit=$(cat "$rea.hit.blast" | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff_tmp -v covcutoff=$covcutoff '{if ($2>=identcutoff && $4>=bitcutoff && $5>=covcutoff) print $0}' | sort -rgk 2,2 | cut -f4)
                        
                        bidihit=$(echo "$back_hit" | grep -Eo "$forward_hit" | sort | uniq | tr '\n' '|' | sed 's/|$//g')
                        #echo forward: $forward_hit
                        #echo bidihit: $bidihit
                        if [ -n "$bidihit" ]; then
                            [[ verbose -ge 1 ]] && echo -e "\t\t\t--> BIDIRECTIONAL hit found <--"
                            grep -E $bidihit $dir/dat/seq/uniprot_sprot.fasta | head -3 | sed -e 's/^/\t\t\t/'
                            is_bidihit=true
                        else
                            is_bidihit=false
                        fi
                    fi
                    
                    if [ -n "$dbhit" ]; then
                        [[ verbose -ge 1 ]] && echo -e '\t\t'Candidate reaction for import: `echo "$dbhit" | wc -w`
                        pwyCand="$pwyCand$dbhit " # remember candidate reaction
                        ((countdb++))
                    else
                        [[ verbose -ge 1 ]] && echo -e '\t\t'NO candidate reaction found for import
                    fi
                    echo "$besthit_all" | awk -v rea="$rea" -v reaName="$reaName" -v ec=$ec -v is_bidihit=$is_bidihit -v dbhit="$dbhit" -v pwy="$pwy" '{print rea"\t"reaName"\t"ec"\t"is_bidihit"\t"$0"\t"pwy"\t""good_blast""\t""NA""\t"dbhit}' >> reactions.tbl
                    ((countex++))
                    countexList="$countexList$rea "
                else
                    someIdentity=$(cat $out | sort -rgk 4,4 | head -1 | cut -f2)
                    someBitscore=$(cat $out | sort -rgk 4,4 | head -1 | cut -f4)
                    someCoverage=$(cat $out | sort -rgk 4,4 | head -1 | cut -f5)
                    somehit_all=$( cat $out | sort -rgk 4,4 | head -1)
                    echo -e "$rea\t$reaName\t$ec\tNA\t$somehit_all\t$pwy\tbad_blast\tNA\t$dbhit" >> reactions.tbl 
                    if [ $subunit_fraction -gt $subunit_cutoff ] || [ $iterations -eq 1 ] ; then
                        [[ verbose -ge 1 ]] && echo -e '\t\t'NO good blast hit"\n\t\t\t(best one: bit=$someBitscore id=$someIdentity cov=$someCoverage)"
                    else
                        [[ verbose -ge 1 ]] && echo -e '\t\t'NO hit because of missing subunits
                    fi 
                    if [[ -n "$dbhit" ]];then
                        pwyNoHitFound="$pwyNoHitFound$dbhit "
                    fi
                fi
            else
                [[ verbose -ge 1 ]] && echo -e '\t\t'NO blast hit
                #echo -e "\t\t$query" 
                if [[ -n "$dbhit" ]];then
                    pwyNoHitFound="$pwyNoHitFound$dbhit "
                    echo -e "$rea\t$reaName\t$ec\tNA\t\t\t\t\t\t\t\t\t$pwy\tno_blast\tNA\t$dbhit" >> reactions.tbl 
                fi
            fi
        else
            [[ verbose -ge 1 ]] && echo -e "\t\tNO sequence data found"
            #echo -e "\t\t$query" 
            #echo -e "\t\t$(basename $query)" 
            ((vague++))
            [[ -n "$dbhit" ]] && pwyVage="$pwyVage$dbhit "
            [[ $keyRea = *"$rea"* ]] && vagueKeyReaFound="$vagueKeyReaFound $rea"
            echo -e "$rea\t$reaName\t$ec\tNA\t\t\t\t\t\t\t\t\t$pwy\tno_seq_data\tNA\t$dbhit" >> reactions.tbl 
        fi
        [[ verbose -ge 2 ]] && echo -e "\t\tCandidate reactions: $dbhit"
    done # pathway
    
    if [ $count -eq 0 ]; then
        completeness=0
    else
        check_vague=$(echo "$vague < $count*$vagueCutoff" | bc) # vague reactions shouldn't make more than certain amount of total reactions
        if [ "$strictCandidates" = false ] && [ $check_vague -eq 1 ] ; then # if vague reaction are considered they should not influence the completeness treshold
            completeness=$(echo "scale=0; 100*($countex+$vague)/$count" | bc)
        else
            completeness=$(echo "scale=0; 100*($countex)/$count" | bc)
        fi
    fi
    if [ $vague -eq 0 ] || [ "$strictCandidates" = true ] || [ $check_vague -eq 0 ]; then
        [[ verbose -ge 1 ]] && echo "Pathway completeness: $countex/$count ($completeness%)"
    else
        [[ verbose -ge 1 ]] && echo "Pathway completeness: ($countex+$vague)/$count ($completeness%) with $vague reactions of unclear state"
    fi

    [[ verbose -ge 1 ]] && echo -e Hits with candidate reactions in database: $countdb/$count
    if [ -n "$keyReaFound" ]; then
        CountKeyReaFound=$(echo $keyReaFound | tr ' ' '\n' |  sort | uniq | wc -l)
    else
        CountKeyReaFound=0
    fi
    
    CountTotalKeyRea=$(echo $keyRea | wc -w)
    CountTotalVagueKeyRea=$(echo $vagueKeyReaFound | wc -w)
    if [[ "$strictCandidates" = false ]] && [[ $CountTotalVagueKeyRea -gt 0 ]]; then
        [[ verbose -ge 1 ]] && echo Key reactions: "$CountKeyReaFound/($CountTotalKeyRea-$CountTotalVagueKeyRea) with $CountTotalVagueKeyRea key reactions of unclear state"
        CountTotalKeyRea=$(echo $CountTotalKeyRea - $CountTotalVagueKeyRea | bc )
    else
        [[ verbose -ge 1 ]] && echo Key reactions: $CountKeyReaFound/$CountTotalKeyRea
    fi
    
    prediction=false
    # add reactions of pathways (even if no blast hit) if above treshold (and no key enzyme is missed)
    cand="$cand$pwyCand " # add all reactions with direct sequence-based evidence
    if [[ $CountTotalKeyRea -gt 0 ]]; then 
        KeyReaFracAvail=$(echo "scale=2; $CountKeyReaFound / $CountTotalKeyRea > 0.5" | bc) # how many key enzymes were found?
    else
        KeyReaFracAvail=1 # no key enzymes are present anyway
    fi

    
    # A) Consider as complete pathway because all reactions are present
    if [[ $completeness -eq 100 ]]; then
        prediction=true
        cand="$cand$pwyVage "
        bestPwy="$bestPwy$name\n"
        pwy_status="full"
    # B) Consider as complete pathway because of completeness treshold (key enzymes should be present too)
    elif [[ $completeness -ge $completenessCutoffNoHints ]] && [[ "$KeyReaFracAvail" -eq 1 ]] && [[ "$strictCandidates" = false ]]; then
        pwy_status="treshold"
        [[ verbose -ge 1 ]] && echo "Consider pathway to be present because of completeness treshold!"
        prediction=true
        cand="$cand$pwyVage$pwyNoHitFound "
        bestPwy="$bestPwy$name ($completeness% completeness, added because of completeness treshhold)\n"
    # C) Consider as complete pathway because of key enzymes (lower treshold)
    elif [[ $CountKeyReaFound -ge 1 ]] && [[ $CountKeyReaFound -eq $CountTotalKeyRea ]] && [[ $completeness -ge $completenessCutoff ]] && [[ "$strictCandidates" = false ]]; then
        pwy_status="keyenzyme"
        [[ verbose -ge 1 ]] && echo "Consider pathway to be present because of key enzyme!"
        prediction=true
        cand="$cand$pwyVage$pwyNoHitFound "
        bestPwy="$bestPwy$name ($completeness% completeness, added because of key enzyme)\n"
    fi
    
    if [[ "$prediction" = true ]];then # update output for reactions when pathway is completete 
        awk -i inplace -v pwy="$pwy" -v pwy_status="$pwy_status" 'BEGIN {OFS=FS="\t"} $13==pwy {$15=pwy_status} 1' reactions.tbl
    fi

    echo -e "$pwy\t$name\t$prediction\t$completeness\t$vague\t$CountTotalKeyRea\t$CountKeyReaFound\t$countexList" >> output.tbl # write down some statistics

done

cand="$(echo $cand | tr ' ' '\n' | sort | uniq | tr '\n' ' ')" # remove duplicates
if [[ verbose -ge 2 ]]; then
    echo -e '\n'Total candidate reactions:
    echo $cand
fi

[[ verbose -ge 1 ]] && echo -e '\n'Pathways found:
[[ verbose -ge 1 ]] && echo -e $bestPwy


# export found reactions 
[[ verbose -ge 1 ]] && echo -e Candidate reactions found: $(echo "$cand" | wc -w) '\n'
echo $cand > newReactions.lst
#cp newReactions.lst $curdir/${fastaID}-$output_suffix-Reactions.lst # not needed anymore
cp output.tbl $curdir/${fastaID}-$output_suffix-Pathways.tbl
[[ -s reactions.tbl ]] && echo "rxn name ec bihit $blast_format pathway status pathway.status dbhit" | tr ' ' '\t' | cat - reactions.tbl | awk '!a[$0]++' > $curdir/${fastaID}-$output_suffix-Reactions.tbl # add header and remove duplicates


# cleaning
[[ -s $tmp_fasta ]] && rm $tmp_fasta


ps -p $$ -o %cpu,%mem,cmd
end_time=`date +%s`
echo Running time: `expr $end_time - $start_time` s.
