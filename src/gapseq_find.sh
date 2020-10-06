#!/bin/bash

start_time=`date +%s`

pathways=""
database="seed"
pwyDatabase="metacyc,custom"
verbose=1
taxonomy="Bacteria"
taxRange="" # taxonomic range for pawthways
bitcutoff=200 # cutoff blast: min bit score
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
use_parallel=true
exhaustive=false
seqSrc=2
anno_genome_cov=false
use_gene_seq=true
stop_on_files_exist=false

usage()
{
    echo "Usage"
    echo "$0 -p keyword / -e ec [-d database] [-t taxonomy] file.fasta."
    echo "  -p keywords such as pathways or subsystems (for example amino,nucl,cofactor,carbo,polyamine)"
    echo "  -e Search by ec numbers (comma separated)"
    echo "  -r Search by enzyme name (colon separated)"
    echo "  -d Database: vmh or seed (default: $database)"
    echo "  -t Taxonomic range for sequences to be downloaded (default: $taxonomy)"
    echo "  -b Bit score cutoff for local alignment (default: $bitcutoff)"
    echo "  -i Identity cutoff for local alignment (default: $identcutoff)"
    echo "  -c Coverage cutoff for local alignment (default: $covcutoff)"
    echo "  -s Strict candidate reaction handling (do _not_ use pathway completeness, key kenzymes and operon structure to infere if imcomplete pathway could be still present (default: $strictCandidates)"
    echo "  -u Suffix used for output files (default: pathway keyword)"
    echo "  -a blast hits back against uniprot enzyme database"
    echo "  -n Consider superpathways of metacyc database"
    echo "  -l Select the pathway database (MetaCyc, KEGG, SEED, all; default: $pwyDatabase)"
    echo "  -o Only list pathways found for keyword; default $onlyList)"
    echo "  -x Do not blast only list pathways, reactions and check for available sequences; default $skipBlast"
    echo "  -q Include sequences of hits in log files; default $includeSeq"

    echo "  -v Verbose level, 0 for nothing, 1 for pathway infos, 2 for full (default $verbose)"
    echo "  -k Do not use parallel"
    echo "  -g Exhaustive search, continue blast even when cutoff is reached (default $exhaustive)"
    echo "  -z Quality of sequences for homology search: 1:only reviewed (swissprot), 2:unreviewed only if reviewed not available, 3:reviewed+unreviewed, 4:only unreviewed (default $seqSrc)"
    echo "  -m Limit pathways to taxonomic range (default $taxRange)"
    echo "  -w Use additional sequences derived from gene names (default $use_gene_seq)"
    echo "  -y Print annotation genome coverage (default $anno_genome_cov)"
    echo "  -j Quit if output files already exist (default $stop_on_files_exist)"

exit 1
}


# paths and eariables
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
uniprotIdentity=0.9 # clustered uniprot database (0.5 or 0.9)
metaPwy=$dir/../dat/meta_pwy.tbl
keggPwy=$dir/../dat/kegg_pwy.tbl
seedPwy=$dir/../dat/seed_pwy.tbl
customPwy=$dir/../dat/custom_pwy.tbl
metaRea=$dir/../dat/meta_rea.tbl
reaDB1=$dir/../dat/vmh_reactions.csv
reaDB2=$dir/../dat/bigg_reactions.tbl
reaDB3=$dir/../dat/seed_reactions_corrected.tsv
reaDB4=$dir/../dat/mnxref_seed.tsv
reaDB5=$dir/../dat/mnxref_seed-other.tsv
reaDB6=$dir/../dat/mnxref_bigg-other.tsv
brenda=$dir/../dat/brenda_ec_edited.csv
seedEC=$dir/../dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv
seedEnzymesNames=$dir/../dat/seed_Enzyme_Name_Reactions_Aliases.tsv
altecdb=$dir/../dat/altec.csv
metaGenes=$dir/../dat/meta_genes.csv


# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "h?p:e:r:d:i:b:c:v:st:nou:al:oxqkgz:m:ywj" opt; do
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
        includeSeq=true
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
    k)
        use_parallel=false
        ;;
    g)
        exhaustive=true
        ;;
    z)
        seqSrc=$OPTARG
        ;;
    m)
        taxRange=$OPTARG
        ;;
    y)
        anno_genome_cov=true
        ;;
    w)
        use_gene_seq=true
        ;;
    j)
        stop_on_files_exist=true
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
        pwyKey="Pathways|seed|kegg"
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

# determine taxonomy
if [ "$taxonomy" == "auto" ]; then
    pred_biom=$($dir/predict_biomass_from16S.sh $fasta)
    if [ "$pred_biom" == "Gram_neg" ] || [ "$pred_biom" == "Gram_pos" ]; then
        taxonomy=Bacteria
    elif [ "$pred_biom" == "Archaea" ]; then
        taxonomy=Archaea
    else
        echo Taxonomy could be predicted automatically.
        exit 1
    fi
    echo Predicted taxonomy: $taxonomy
fi

# squence directory
export LC_NUMERIC="en_US.UTF-8"
seqpath=$dir/../dat/seq/$taxonomy
seqpath_user=$dir/../dat/seq/$taxonomy/user
mkdir -p $seqpath/rev $seqpath/unrev $seqpath_user

# download sequences if needed
if [[ ! -f $seqpath/rev/sequences.tar.gz  ]] || [[ ! -f $seqpath/unrev/sequences.tar.gz ]] || [[ ! -f $seqpath/rxn/sequences.tar.gz ]]; then
    $dir/update_sequences.sh $taxonomy
fi


if [ -n "$ecnumber" ] || [ -n "$reaname" ]; then
    # create dummpy pwy template for given ec number
    if [[ -z "$ecnumber" ]]; then
        rea_count=$(echo $reaname | tr ';' '\n' | wc -l)
        ecnumber=$(echo $reaname | grep -o ";" | tr ';' ',') # get dummy empty comma seperated ec numbers
    elif [[ -z "$reaname" ]]; then
        rea_count=$(echo $ecnumber | tr ',' '\n' | wc -l)
        reaname=$(echo $ecnumber | grep -o "," | tr -d '\n' | tr ',' ';') # get dummy empty colon seperated reaction names
    else  
        rea_count=$(echo $ecnumber | tr ',' '\n' | wc -l)
    fi
    rea_id=$(seq 1 $rea_count | awk '{print "reaction"$1}' | tr '\n' ',' | sed 's/,$//g')
    pwyKey="custom"
    #pwyDB=$(echo -e "custom\t$ecnumber\t\t\t\t${rea_id::-1}\t$ecnumber\t\t$reaname") # slicing is incompatible?
    pwyDB=$(echo -e "custom\t$ecnumber\t\t\t\t${rea_id}\t$ecnumber\t\t$reaname")
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
    dupli=$(cat allPwy | cut -f1 | sort | uniq -d | tr -d 'id' | sed '/^$/d')
    if [ -n "$dupli" ]; then
        [[ verbose -ge 1 ]] && echo Duplicated pathway IDs found: $dupli will only use $customPwy
        dupli_search=$(echo "$dupli" | sed 's/|/\\|/g' |tr '\n' '|' | rev | cut -c2- | rev)
        [[ verbose -ge 1 ]] && echo "$dupli_search"
        cat allPwy | grep -wEv "$dupli_search" > allPwy.tmp
        cat $customPwy | grep -wE "$dupli_search" >> allPwy.tmp
        mv allPwy.tmp allPwy
        #cat allPwy | grep -wE "$dupli_search"
    fi
    pwyDB=$(cat allPwy | grep -wEi $pwyKey)
    [[ "$noSuperpathways" = true ]] && pwyDB=$(echo "$pwyDB" | grep -v 'Super-Pathways')
    [ -z "$ecnumber" ] && [ -z "$pwyDB" ] && { echo "No pathways found for key $pwyKey"; exit 1; }
fi
[ -z "$output_suffix" ] && output_suffix=$pathways


# check to quit if output files already exist
if [[ "$stop_on_files_exist" = true ]] && [[ -f $curdir/${fastaID}-$output_suffix-Pathways.tbl ]] && [[ -f $curdir/${fastaID}-$output_suffix-Reactions.tbl ]]; then
    echo Pathway and reaction output files already exist, exiting...
    exit 0
fi


# function to get database hits for ec number
getDBhit(){
    kegg=$(grep -wFe "|$rea" $metaRea | awk -F "\t" {'print $5'})

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
    if [ "$database" == "vmh" ]; then
        dbhit="$dbhit $(grep -wFe "$rea" $reaDB2 | awk '{print $1}')"
    fi

    # 5) match reaction using mnxref namespace
    if [ "$database" == "seed" ]; then
        dbhit="$dbhit $(grep -wFe "|$rea" $reaDB5 | awk '{print $2}')"
    elif [ "$database" == "vmh" ]; then
        dbhit="$dbhit $(grep -wFe "|$rea" $reaDB6 | awk '{print $2}')"
    fi

    # 6) match reaction using custom enzyme-name - seedID mapping
    if [ "$database" == "seed" ] & [ "$reaName" != "" ]; then
        dbhit="$dbhit $(grep -wFe "$reaName" $seedEnzymesNames | awk -F '\t' ' {print $1}')"
    fi

    [ "$dbhit" == " " ] && dbhit=""
}





# create blast database
makeblastdb -in $fasta -dbtype nucl -out orgdb >/dev/null


cand=""     #list of candidate reactions to be added
bestPwy=""  # list of found pathways
echo -e "ID\tName\tPrediction\tCompleteness\tVagueReactions\tKeyReactions\tKeyReactionsFound\tReactionsFound" > output.tbl # pahtway statistics file

#taxRange=Proteobacteria
if [ -n "$taxRange" ]; then
    validTax=$(grep -i $taxRange $dir/../dat/taxonomy.tbl | cut -f1 | tr '\n' '|' | sed 's/.$//')
    pwyDB_new=$(echo "$pwyDB" | grep -wE `echo TAX-$validTax`)
    pwyDB_old=$(echo "$pwyDB" | awk -F '\t' 'BEGIN {OFS=FS="\t"} $5=="" {print $0}')
    pwyDB=$(echo "$pwyDB_new""$pwyDB_old")
    [[ -z "$pwyDB" ]] && { echo "No pathways found"; exit 0; }
fi

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
    spontRea=$(echo "$line" | awk -F "\t" '{print $14}' | tr ',' ' ')
    pwyHierarchy=$(echo "$line" | awk -F "\t" '{print $4}' | sed 's/|\|THINGS\|Generalized-Reactions\|Pathways\|FRAMES//g' | sed 's/,,//g')
    reaNr=$(echo $ecs | tr "," "\n" | wc -l)
    [[ verbose -ge 1 ]] && echo -e '\n'$i/$pwyNr: Checking for pathway $pwy $name with $reaNr reactions
    [[ verbose -ge 1 ]] && echo "($pwyHierarchy)"
    [[ "$onlyList" = true ]] && { continue; }
    for j in `seq 1 $reaNr`
    #for ec in $(echo $ecs | tr "," "\n")
    do 
        dbhit=""
        ec=$(echo $ecs | awk -v j=$j -F ',' '{print $j}')
        rea=$(echo $reaids | awk -v j=$j -F ',' '{print $j}')
        reaName=$(echo $reaNames | awk -v j=$j -F ';' '{print $j}' | tr -d '|')
        re="([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)"
        EC_test=$(if [[ $ec =~ $re ]]; then echo ${BASH_REMATCH[1]}; fi) # check if not trunked ec number (=> too many hits)
        geneName=$(cat $metaGenes | awk -v rea=$rea -v pwy=$pwy -F ',' '$1~rea && $3==pwy {print $2}')
        geneRef=$(cat $metaGenes | awk -v rea=$rea -v pwy=$pwy -F ',' '$1~rea && $3==pwy {print $5}')
        [[ verbose -ge 1 ]] && echo -e "\t$j) $rea $reaName $ec" $geneName
        [[ -z "$rea" ]] && { continue; }
        [[ -n "$ec" ]] && [[ -n "$reaName" ]] && [[ -n "$EC_test" ]] && { is_exception=$(grep -Fw -e "$ec" -e "$reaName" $dir/../dat/exception.tbl | wc -l); }
        ( [[ -z "$ec" ]] || [[ -z "$EC_test" ]] ) && [[ -n "$reaName" ]] && { is_exception=$(grep -Fw "$reaName" $dir/../dat/exception.tbl | wc -l); }
        [[ -n "$ec" ]] && [[ -z "$reaName" ]] && [[ -n "$EC_test" ]] && { is_exception=$(grep -Fw "$ec" $dir/../dat/exception.tbl | wc -l); }
        if [[ $is_exception -gt 0 ]] && [[ $identcutoff -lt $identcutoff_exception ]];then # take care of similair enzymes with different function
            identcutoff_tmp=$identcutoff_exception
            [[ verbose -ge 1 ]] && echo -e "\t\tUsing higher identity cutoff for $rea"
        else
            identcutoff_tmp=$identcutoff
        fi
        getDBhit # get db hits for this reactions
        dbhit="$(echo $dbhit | tr ' ' '\n' | sort | uniq | tr '\n' ' ')" # remove duplicates
        #dbhit=$($dir/getDBhit.sh "$rea" "$reaName" "$ec" "$database" "$EC_test")
        if [[ $spontRea = *"$rea"* ]]; then # detect spontaneous reactions
            [[ verbose -ge 1 ]] && echo -e '\t\t--> Spontaneous reaction <--'
            if [ "$includeSeq" = true ]; then
                echo -e "$rea\t$reaName\t$ec\tNA\t\t\t\t\t\t\t\t\t\t$pwy\tspontaneous\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
            else
                echo -e "$rea\t$reaName\t$ec\tNA\t\t\t\t\t\t\t\t\t$pwy\tspontaneous\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
            fi
            continue
        fi
        ((count++))
        if [[ -n "$EC_test" ]]; then
            # check if sequence is not available => try to download
            if [ ! -f $seqpath/rev/$ec.fasta ]; then
                [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading reviewed sequences for: $ec 
                $dir/uniprot.sh -e "$ec" -t "$taxonomy" -i $uniprotIdentity >/dev/null
            fi
            if [ ! -f $seqpath/unrev/$ec.fasta ] && [ $seqSrc -gt 1 ]; then 
                [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading unreviewed sequences for: $ec 
                 $dir/uniprot.sh -u -e "$ec" -t "$taxonomy" -i $uniprotIdentity >/dev/null
            fi
            if [ -s "$seqpath_user/$ec.fasta" ]; then
                [[ verbose -ge 1 ]] && echo -e "\t\t--> Found user defined sequence file <--"
                query=$seqpath_user/$ec.fasta
            elif [ $seqSrc -eq 1 ]; then
                query=$seqpath/rev/$ec.fasta
            elif [ $seqSrc -eq 2 ]; then
                query=$seqpath/rev/$ec.fasta
                [[ ! -s $query ]] && query=$seqpath/unrev/$ec.fasta
            elif [ $seqSrc -eq 3 ]; then
                query=$(mktemp)
                cat $seqpath/rev/$ec.fasta $seqpath/unrev/$ec.fasta > $query
            elif [ $seqSrc -eq 4 ]; then
                query=$seqpath/unrev/$ec.fasta
            fi
            
            # if an alternative ec numbers exists and has additional sequence data merge both files
            EC_test2=$(if [[ $altec =~ $re ]]; then echo ${BASH_REMATCH[1]}; fi) # check if not trunked ec number (=> too many hits)
            if [[ -n "$EC_test2" ]]; then
                [[ verbose -ge 1 ]] && echo -e "\t\t--> Found alternative EC number:" $altec "<--"
                # check if sequence is not available => try to download
                #[[ `echo "$altec" | wc -l` -gt 1 ]] && exit 1 # multiple alternative EC
                query_alt_all=$(mktemp)
                for aec in $altec; do
                    if [ ! -f $seqpath/rev/$aec.fasta ]; then
                        [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading reviewed sequences for: $aec 
                        $dir/uniprot.sh -e "$aec" -t "$taxonomy" -i $uniprotIdentity >/dev/null
                    fi
                    if [ ! -f $seqpath/unrev/$aec.fasta ] && [ $seqSrc -gt 1 ]; then 
                        [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading unreviewed sequences for: $aec 
                        $dir/uniprot.sh -u -e "$aec" -t "$taxonomy" -i $uniprotIdentity >/dev/null
                    fi
                    if [ -s "$seqpath_user/$aec.fasta" ]; then
                        [[ verbose -ge 1 ]] && echo -e "\t\t--> Found user defined sequence file <--"
                        query_alt=$seqpath_user/$aec.fasta
                    elif [ $seqSrc -eq 1 ]; then
                        query_alt=$seqpath/rev/$aec.fasta
                    elif [ $seqSrc -eq 2 ]; then
                        query_alt=$seqpath/rev/$aec.fasta
                        [[ ! -s $query_alt ]] && query_alt=$seqpath/unrev/$aec.fasta
                    elif [ $seqSrc -eq 3 ]; then
                        query_alt=$(mktemp)
                        cat $seqpath/rev/$aec.fasta $seqpath/unrev/$aec.fasta > $query_alt
                    elif [ $seqSrc -eq 4 ]; then
                        query_alt=$seqpath/unrev/$aec.fasta
                    fi
                    cat $query_alt >> $query_alt_all
                done
                #merge sequence data
                if [[ -s $query_alt_all ]]; then
                    [[ verbose -ge 1 ]] && { echo -e "\t\tMerge sequence data from `basename $query` and" $altec;  }
                    query_merge=$(mktemp)
                    cat $query $query_alt_all | awk '/^>/{f=!d[$1];d[$1]=1}f' > $query_merge # no duplicates
                    query=$query_merge
                fi
            fi
        fi
        
        # if no EC number is available or no sequence was found for EC number then use reaction name instead for sequence search
        if [[ -n "$reaName" ]] && ( [[ -z "$EC_test" ]] || [[ ! -s "$query" ]] );then
            reaNameHash=$(echo -n "$reaName" | md5sum | awk '{print $1}')
            # check if sequence is not available => try to download
            if [ ! -f $seqpath/rev/$reaNameHash.fasta ]; then
                [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading reviewed sequences for: $reaName "\n\t\t(hash: $reaNameHash)" 
                $dir/uniprot.sh -r "$reaName" -t "$taxonomy" -i $uniprotIdentity >/dev/null
            fi
            if [ ! -f $seqpath/unrev/$reaNameHash.fasta ] && [ $seqSrc -gt 1 ]; then 
                [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading unriewed sequences for: $reaName "\n\t\t(hash: $reaNameHash)" 
                $dir/uniprot.sh -u -r "$reaName" -t "$taxonomy" -i $uniprotIdentity >/dev/null
            fi
            if [ -s "$seqpath_user/$reaNameHash.fasta" ]; then
                [[ verbose -ge 1 ]] && echo -e "\t\t--> Found user defined sequence file <--"
                query=$seqpath_user/$reaNameHash.fasta
            elif [ $seqSrc -eq 1 ]; then
                query=$seqpath/rev/$reaNameHash.fasta
            elif [ $seqSrc -eq 2 ]; then
                query=$seqpath/rev/$reaNameHash.fasta
                [[ ! -s $query ]] && query=$seqpath/unrev/$reaNameHash.fasta
            elif [ $seqSrc -eq 3 ]; then
                query=$(mktemp)
                cat $seqpath/rev/$reaNameHash.fasta $seqpath/unrev/$reaNameHash.fasta > $query
            elif [ $seqSrc -eq 4 ]; then
                query=$seqpath/unrev/$reaNameHash.fasta
            fi
        fi

        # sequence by gene name
        if [[ -n "$geneName" ]] && [[ -n "$geneRef" ]] && [[ "$use_gene_seq" = true ]]; then
            if [ ! -f $seqpath/rxn/$rea.fasta ]; then
                [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading sequences for: $geneRef
                reaSeqTmp=$(mktemp)
                for gr in $geneRef
                do
                    $dir/uniprot.sh -d $gr -t "$taxonomy" -i $uniprotIdentity >/dev/null
                    if [ -f $seqpath/rxn/$gr.fasta ]; then
                        cat $seqpath/rxn/$gr.fasta >> $reaSeqTmp
                        rm $seqpath/rxn/$gr.fasta
                    fi
                done
                if [ -s "$reaSeqTmp" ]; then
                   mv $reaSeqTmp $seqpath/rxn/$rea.fasta
                else
                    touch $seqpath/rxn/$rea.fasta # create empty file if no gene seq data is found to avoid reoccuring download attempt
                fi
            fi
            
            if [ -s "$seqpath_user/$rea.fasta" ]; then
                [[ verbose -ge 1 ]] && echo -e "\t\t--> Found user defined sequence file <--"
                query_gene=$seqpath_user/$rea.fasta
            else
                query_gene=$seqpath/rxn/$rea.fasta
            fi
            #merge sequence data
            if [[ -s $query_gene ]]; then
                [[ verbose -ge 1 ]] && { echo -e "\t\tMerge sequence data from `basename $query` and" `basename $query_gene`;  }
                query_merge2=$(mktemp)
                cat $query $query_gene | awk '/^>/{f=!d[$1];d[$1]=1}f' > $query_merge2 # no duplicates
                query=$query_merge2
            fi
        
        fi
        
        [[ "$skipBlast" = true ]] && { continue; }

        if [ -s $query ]; then
            [[ verbose -ge 1 ]] && echo -e "\t\t$query (`cat $query | grep ">" | wc -l` sequences)"
            #out=$rea.blast #$ec.blast
            query_id=$(basename $query)
            out="${query_id%.fasta}".blast
            subunits_found=0
            subunits_undefined_found=0
            subunit_prescan=0
            iteractions=0
            subunits_former_run=false
            if [ ! -f $out ]; then # check if there is a former hit
                subunit_prescan=$(cat $query | sed -n 's/^>//p' | grep -E 'subunit|chain|polypeptide|component' | wc -l) # prescan if subunits can be found because detection script is time intensive
                if [ $subunit_prescan -gt 0 ]; then
                    Rscript $dir/complex_detection.R $query subunit_tmp.fasta # set new fasta header with consistent subunit classification (avoid mix of arabic,latin and greek numbers)
                    query=$(readlink -f subunit_tmp.fasta)
                fi

                #subunits=$(cat $query | sed -n 's/^>//p' | grep -oE 'subunit [0-9]|(alpha|beta|gamma|delta|epsilon) subunit' | sort | uniq) # check for subunits
                subunits=$(cat $query | sed -n 's/^>//p' | grep -oE 'Subunit \w+$' | sort | uniq) # check for subunits
                subunits_count=$(echo -e "$subunits"| wc -l) 
                undefined=$(cat $query | sed -n 's/^>//p' | grep -Ev 'Subunit \w+$' | sort | uniq) # check for sequences which do not follow regular expression => will be treated as other (i.e. one additional subunit)
                [[ -n "$undefined" ]] && subunits=$(echo -e "$subunits\nSubunit undefined" | sed '/^$/d') # add default case for undefined subunits
                [[ verbose -ge 1 ]] && [[ $subunits_count -gt 1 ]] && echo -e '\t\t'check subunits: $subunits_count
                iterations=$(echo -e "$subunits"| wc -l) # every subunit will get a own iteration
                for iter in `seq $iterations`
                do
                    if [ $iterations -gt 1 ]; then 
                        # apt install exonerate
                        fastaindex $query query.idx
                        subunit_id=$(echo "$subunits" | sed -n ${iter}p)
                        #echo $subunit_id
                        if [ "$subunit_id" == "Subunit undefined" ];then
                            subunit_id2=$(echo "$subunits" | tr '\n' '|' | sed 's/|$//g') # inverse search
                            cat $query | sed -n 's/^>//p' | grep -Ev "$subunit_id2" | awk '{print $1}' | sed 's/^>//g' > query_subunit_header
                        else
                            cat $query | grep "$subunit_id" | awk '{print $1}' | sed 's/^>//g' > query_subunit_header
                        fi
                        #echo -e $iter "\n"
                        fastafetch -f $query -i query.idx -Fq <(sort -u query_subunit_header) > query_subunit.fasta
                        rm query.idx query_subunit_header
                    else
                        cp $query query_subunit.fasta
                    fi
                    #csplit -s -z query_subunit.fasta '/>/' '{*}' # split multiple fasta file and skip further testing if a significant hit is found
                    $dir/fasta-splitter.pl --n-parts 10 query_subunit.fasta >/dev/null
                    #for q in `ls xx*`
                    subunits_found_old=$subunits_found
                    for q in `ls query_subunit.part-*.fasta`
                    do
                        if ! [ -x "$(command -v parallel)" ] || [ "$use_parallel" = false ]; then # try to use parallelized version
                            tblastn -db orgdb -query $q -qcov_hsp_perc $covcutoff -outfmt "6 $blast_format" > query.blast
                        else
                            cat $q | parallel --gnu --will-cite --block 50k --recstart '>' --pipe tblastn -db orgdb -qcov_hsp_perc $covcutoff -outfmt \'"6 $blast_format"\' -query - > query.blast
                        fi
                        cat query.blast >> $out
                        bhit=$(cat query.blast | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff_tmp -v covcutoff=$covcutoff '{if ($2>=identcutoff && $4>=bitcutoff && $5>=covcutoff) print $0}')
                        if [ -n "$bhit" ]; then
                            bestsubunithit=$(echo "$bhit" | sort -rgk 4,4)
                            [[ verbose -ge 1 ]] && [[ $iterations -gt 1 ]] && echo "$bestsubunithit" | head -1 | cut -f1 | grep -f - $query | sed "s/^/\t\t\t$subunit_id hit: /" 
                            ((subunits_found++))
                            [[ "$subunit_id" == "Subunit undefined" ]] && ((subunits_undefined_found++))
                            [[ $iterations -gt 1 ]] && echo "$bestsubunithit" | awk -v exception="$is_exception" -v subunit="$subunit_id" -v rea="$rea" -v reaName="$reaName" -v ec=$ec -v dbhit="$dbhit" -v pwy="$pwy" '{print rea"\t"reaName"\t"ec"\t"NA"\t"$0"\t"pwy"\t""good_blast""\t""NA""\t"dbhit"\t"subunit"\t"exception"\t""NA"}' >> reactions.tbl
                            [[ "$exhaustive" = false ]] && break 
                        fi
                        rm query.blast
                    done
                    if [[ $subunits_found -eq $subunits_found_old ]] && [[ $iterations -gt 1 ]]; then
                        if [ "$includeSeq" = true ]; then
                            echo -e "$rea\t$reaName\t$ec\tNA\t\t\t\t\t\t\t\t\t\t$pwy\tno_blast\tNA\t$dbhit\t$subunit_id\t$is_exception\tNA" >> reactions.tbl # subunit not found 
                        else
                            echo -e "$rea\t$reaName\t$ec\tNA\t\t\t\t\t\t\t\t\t$pwy\tno_blast\tNA\t$dbhit\t$subunit_id\t$is_exception\tNA" >> reactions.tbl # subunit not found 
                        fi
                    fi
                    rm query_subunit.part-*.fasta*
                done
                [[ $iterations -gt 1 ]] && [[ verbose -ge 1 ]] &&  echo -e '\t\t'total subunits found: `echo $subunits_found - $subunits_undefined_found | bc` / $subunits_count
                #[[ $iterations -gt 1 ]] && [[ verbose -ge 1 ]] && [[ $subunits_undefined_found -eq 1 ]] && echo -e '\t\tUndefined subunit found' 
                echo -e $out'\t'$subunits_found'\t'$iterations'\t'$subunits_count'\t'$subunits_undefined_found >> subunits.log # save subunits found
            else
                # get subunit fraction from former run
                subunits_former_run=true
                subunits_found=$(cat subunits.log | awk -F "\t" -v out=$out '$1==out { print $2 }')
                iterations=$(cat subunits.log | awk -F "\t" -v out=$out '$1==out { print $3 }')
                subunits_count=$(cat subunits.log | awk -F "\t" -v out=$out '$1==out { print $4 }')
                subunits_undefined_found=$(cat subunits.log | awk -F "\t" -v out=$out '$1==out { print $5 }')
                #echo test:$subunits_found $iterations
            fi
            if [ -s $out ]; then
                bhit=$(cat $out | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff_tmp -v covcutoff=$covcutoff '{if ($2>=identcutoff && $4>=bitcutoff && $5>=covcutoff) print $0}')
                subunit_fraction=$(echo "100*$subunits_found/$subunits_count" | bc)
                [[ $subunit_fraction -eq $subunit_cutoff ]] && [[ $subunits_undefined_found -eq 1 ]] && [[ verbose -ge 1 ]] && echo -e '\t\t\tUndefined subunit caused that threshold is passed' # undefined subunit can have bonus effect
                if [ -n "$bhit" ] && ( [ $subunit_fraction -gt $subunit_cutoff ] || ( [ $subunit_fraction -eq $subunit_cutoff ] && [ $subunits_undefined_found -eq 1 ] ) ); then
                    bestIdentity=$(echo "$bhit" | sort -rgk 4,4 | head -1 | cut -f2)
                    bestBitscore=$(echo "$bhit" | sort -rgk 4,4 | head -1 | cut -f4)
                    bestCoverage=$(echo "$bhit" | sort -rgk 4,4 | head -1 | cut -f5)
                    besthit_all=$(echo "$bhit" | sort -rgk 4,4)
                    bhit_count=$(echo "$bhit" | wc -l)
                    [[ verbose -ge 1 ]] && echo -e '\t\t'Blast hit \(${bhit_count}x\)
                    [[ verbose -ge 1 ]] && [[ $iterations -le 1 ]] && echo "$besthit_all" | head -3 | awk '{print "\t\t\tbit="$4 " id="$2 " cov="$5 " hit="$1}' # only for non-subunit hits
                    # check if key reactions of pathway
                    if [[ $keyRea = *"$rea"* ]]; then
                        [[ verbose -ge 1 ]] && echo -e '\t\t--> KEY reaction found <--'
                        keyReaFound="$keyReaFound $rea"
                    fi
                    [[ $iterations -le 1 ]] && echo "$besthit_all" | awk -v exception="$is_exception" -v rea="$rea" -v reaName="$reaName" -v ec=$ec -v dbhit="$dbhit" -v pwy="$pwy" '{print rea"\t"reaName"\t"ec"\t"NA"\t"$0"\t"pwy"\t""good_blast""\t""NA""\t"dbhit"\t""NA""\t"exception"\t""NA"}' >> reactions.tbl # only for non-subunit hits
                    awk -v rea="$rea" -v status="1" 'BEGIN {OFS=FS="\t"} $1==rea {$19=status} 1' reactions.tbl > reactions.tmp.tbl && mv reactions.tmp.tbl reactions.tbl # change protein complex status for all subunits found
                    
#blast hit back to uniprot enzyme database
                    if [ "$blast_back" = true ]; then
                        [[ verbose -ge 1 ]] && echo -e "\t\tBlast best hits against uniprot db:"
                        for iter in `seq $iterations`
                        do
                            subiter_tmp=$(echo Subunit $iter)
                            #echo $iterations $rea $subiter_tmp
                            if [ $iterations -le 1 ]; then
                                subiter_log=$(cat reactions.tbl | awk -F "\t" -v rea="$rea" '$1==rea { print $0 }')
                            else
                                subiter_log=$(cat reactions.tbl | awk -F "\t" -v rea="$rea" -v subiter_tmp="$subiter_tmp" '$1==rea && $18==subiter_tmp { print $0 }')
                            fi
                            echo "$subiter_log" | awk -F "\t" '{print ">"$5"\n"$13}' > $rea.hit.$iter.fasta
                            blastp -db $dir/../dat/seq/uniprot_sprot -query "$rea.hit.$iter.fasta" -outfmt '6 pident bitscore qcovs sseqid qseqid' > $rea.hit.blast 2>/dev/null

                            #forward_hit=$(echo "$bhit" | sort -rgk 4,4 | cut -f1 | sed 's/UniRef90_//g' | sort | uniq | tr '\n' '|' | sed 's/|$//g')
                            forward_hit=$(echo "$subiter_log" | awk -F "\t" '{print $5}' | sort -rgk 4,4 | cut -f1 | sed 's/UniRef90_//g' | sort | uniq | tr '\n' '|' | sed 's/|$//g')
                            back_hit=$(cat "$rea.hit.blast" | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff_tmp -v covcutoff=$covcutoff '{if ($2>=identcutoff && $4>=bitcutoff && $5>=covcutoff) print $0}' | sort -rgk 2,2 | cut -f4)
                            
                            bidihit=$(echo "$back_hit" | grep -Eo "$forward_hit" | sort | uniq | tr '\n' '|' | sed 's/|$//g')
                            #echo forward: $forward_hit
                            #echo bidihit: $bidihit
                            if [ -n "$bidihit" ]; then
                                if [ $iterations -le 1 ]; then
                                    [[ verbose -ge 1 ]] && echo -e "\t\t\t--> BIDIRECTIONAL hit found <--"
                                else
                                    [[ verbose -ge 1 ]] && echo -e "\t\t\t--> BIDIRECTIONAL hit found for $subiter_tmp <--"
                                fi
                                grep -E $bidihit $dir/../dat/seq/uniprot_sprot.fasta | head -3 | sed -e 's/^/\t\t\t    /'
                                is_bidihit=true
                            else
                                is_bidihit=false
                            fi
                            if [ $iterations -le 1 ]; then
                                awk -v rea="$rea" -v is_bidihit="$is_bidihit" 'BEGIN {OFS=FS="\t"} $1==rea {$4=is_bidihit} 1' reactions.tbl > reactions.tmp.tbl && mv reactions.tmp.tbl reactions.tbl # change bidirectional status
                            else
                                awk -v rea="$rea" -v is_bidihit="$is_bidihit" -v subiter_tmp="$subiter_tmp" 'BEGIN {OFS=FS="\t"} $1==rea && $18==subiter_tmp {$4=is_bidihit} 1' reactions.tbl > reactions.tmp.tbl && mv reactions.tmp.tbl reactions.tbl # change bidirectional status
                            fi
                        done
                    fi
                    
                    if [ -n "$dbhit" ]; then
                        [[ verbose -ge 1 ]] && echo -e '\t\t'Candidate reaction for import: `echo "$dbhit" | wc -w`
                        pwyCand="$pwyCand$dbhit " # remember candidate reaction
                        ((countdb++))
                    else
                        [[ verbose -ge 1 ]] && echo -e '\t\t'NO candidate reaction found for import
                    fi

                    ((countex++))
                    countexList="$countexList$rea "
                else
                    someIdentity=$(cat $out | sort -rgk 4,4 | head -1 | cut -f2)
                    someBitscore=$(cat $out | sort -rgk 4,4 | head -1 | cut -f4)
                    someCoverage=$(cat $out | sort -rgk 4,4 | head -1 | cut -f5)
                    somehit_all=$( cat $out | sort -rgk 4,4 | head -1)
                    if [ $subunit_fraction -gt $subunit_cutoff ] || [ $iterations -eq 1 ] ; then
                        [[ verbose -ge 1 ]] && echo -e '\t\t'NO good blast hit"\n\t\t\t(best one: bit=$someBitscore id=$someIdentity cov=$someCoverage)"
                        echo -e "$rea\t$reaName\t$ec\tNA\t$somehit_all\t$pwy\tbad_blast\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
                    else
                        [[ verbose -ge 1 ]] && echo -e '\t\t'NO hit because of missing subunits
                        if [[ "$subunits_former_run" = true ]];then # log also subunits from former run
                            tmp_log=$(cat reactions.tbl | awk -F '\t' -v rea="$rea" -v reaName="$reaName" -v ec="$ec" -v pwy="$pwy" '{OFS=FS} $1==rea && $2==reaName && $3==ec {$13=pwy; print}')
                            echo "$tmp_log" >> reactions.tbl
                        fi
                    fi 
                    if [[ -n "$dbhit" ]];then
                        pwyNoHitFound="$pwyNoHitFound$dbhit "
                    fi
                fi
            else
                [[ verbose -ge 1 ]] && { echo -e '\t\t'NO blast hit; }
                if [[ -n "$dbhit" ]];then
                    pwyNoHitFound="$pwyNoHitFound$dbhit "
                    if [ "$includeSeq" = true ]; then
                        echo -e "$rea\t$reaName\t$ec\tNA\t\t\t\t\t\t\t\t\t\t$pwy\tno_blast\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
                    else
                        echo -e "$rea\t$reaName\t$ec\tNA\t\t\t\t\t\t\t\t\t$pwy\tno_blast\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
                    fi
                fi
            fi
        else
            [[ verbose -ge 1 ]] && { echo -e "\t\tNO sequence data found"; echo -e "\t\t$query"; }
            ((vague++))
            [[ -n "$dbhit" ]] && pwyVage="$pwyVage$dbhit "
            [[ $keyRea = *"$rea"* ]] && vagueKeyReaFound="$vagueKeyReaFound $rea"
            if [ "$includeSeq" = true ]; then
                echo -e "$rea\t$reaName\t$ec\tNA\t\t\t\t\t\t\t\t\t\t$pwy\tno_seq_data\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
            else
                echo -e "$rea\t$reaName\t$ec\tNA\t\t\t\t\t\t\t\t\t$pwy\tno_seq_data\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
            fi
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
        #KeyReaFracAvail=$(echo "scale=2; $CountKeyReaFound / $CountTotalKeyRea > 0.5" | bc) # how many key enzymes were found?
        KeyReaFracAvail=$(echo "scale=2; $CountKeyReaFound / $CountTotalKeyRea == 1" | bc) # how many key enzymes were found?
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
        #awk -i inplace -v pwy="$pwy" -v pwy_status="$pwy_status" 'BEGIN {OFS=FS="\t"} $13==pwy {$15=pwy_status} 1' reactions.tbl
        awk -v pwy="$pwy" -v pwy_status="$pwy_status" 'BEGIN {OFS=FS="\t"} $13==pwy {$15=pwy_status} 1' reactions.tbl > reactions.tmp.tbl && mv reactions.tmp.tbl reactions.tbl
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
[[ -s reactions.tbl ]] && echo "rxn name ec bihit $blast_format pathway status pathway.status dbhit complex exception complex.status" | tr ' ' '\t' | cat - reactions.tbl | awk '!a[$0]++' > $curdir/${fastaID}-$output_suffix-Reactions.tbl # add header and remove duplicates

# print annotation genome coverage
[[ verbose -ge 1 ]] && [[ "$anno_genome_cov" = true ]] && Rscript $dir/coverage.R $fasta $curdir/${fastaID}-$output_suffix-Reactions.tbl 

# cleaning
[[ -s $tmp_fasta ]] && rm $tmp_fasta


ps -p $$ -o %cpu,%mem,cmd
end_time=`date +%s`
echo Running time: `expr $end_time - $start_time` s.
