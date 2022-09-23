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
update_manually=false
user_temp=false
force_offline=false
input_mode="auto"
n_threads=`grep -c ^processor /proc/cpuinfo`

usage()
{
    echo "Usage"
    echo "$0 -p <keyword> / -e <EC> [-d <database>] [-t <taxonomy>] file.fasta"
    echo "  -p keywords such as pathways or subsystems (for example amino,nucl,cofactor,carbo,polyamine)"
    echo "  -e Search by ec numbers (comma separated)"
    echo "  -r Search by enzyme name (colon separated)"
    echo "  -d Database: vmh or seed (default: $database)"
    echo "  -t Taxonomic range for reference sequences to be used. (Bacteria, Archaea, auto; default: $taxonomy). See Details."
    echo "  -b Bit score cutoff for local alignment (default: $bitcutoff)"
    echo "  -i Identity cutoff for local alignment (default: $identcutoff)"
    echo "  -c Coverage cutoff for local alignment (default: $covcutoff)"
    echo "  -s Strict candidate reaction handling (do _not_ use pathway completeness, key kenzymes and operon structure to infere if imcomplete pathway could be still present (default: $strictCandidates)"
    echo "  -u Suffix used for output files (default: pathway keyword)"
    echo "  -a blast hits back against uniprot enzyme database"
    echo "  -n Consider superpathways of metacyc database"
    echo "  -l Select the pathway database (MetaCyc, KEGG, SEED, all; default: $pwyDatabase)"
    echo "  -o Only list pathways found for keyword (default: $onlyList)"
    echo "  -x Do not blast only list pathways, reactions and check for available sequences (default: $skipBlast)"
    echo "  -q Include sequences of hits in log files (default: $includeSeq)"
    echo "  -v Verbose level, 0 for nothing, 1 for pathway infos, 2 for full (default: $verbose)"
    echo "  -k Do not use parallel (Deprecated: use '-K 1' instead to disable multi-threading.)"
    echo "  -g Exhaustive search, continue blast even when cutoff is reached (default: $exhaustive)"
    echo "  -z Quality of sequences for homology search: 1:only reviewed (swissprot), 2:unreviewed only if reviewed not available, 3:reviewed+unreviewed, 4:only unreviewed (default: $seqSrc)"
    echo "  -m Limit pathways to taxonomic range (default: $taxRange)"
    echo "  -w Use additional sequences derived from gene names (default: $use_gene_seq)"
    echo "  -y Print annotation genome coverage (default: $anno_genome_cov)"
    echo "  -j Quit if output files already exist (default: $stop_on_files_exist)"
    echo "  -U Do not use gapseq sequence archive and update sequences from uniprot manually (very slow) (default: $update_manually)"
    echo "  -T Set user-defined temporary folder (default: $user_temp)"
    echo "  -O For offline mode (default: $force_offline)"
    echo "  -M Input genome mode. Either 'nucl', 'prot', or 'auto' (default $input_mode)"
    echo "  -K Number of threads for sequence alignments. If option is not provided, number of available CPUs will be automatically determined."
    echo ""
    echo "Details:"
    echo "\"-t\": if 'auto', gapseq tries to predict if the organism is Bacteria or Archaea based on the provided genome sequence. The prediction is based on the 16S rRNA gene sequence using a classifier that was trained on 16S rRNA genes from organisms with known Gram-staining phenotype. In case no 16S rRNA gene was found, a k-mer based classifier is used instead."

exit 1
}


# paths and eariables
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
script_name=$(basename -- "$0")
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

function join_by { local IFS="$1"; shift; echo "$*"; }

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "h?p:e:r:d:i:b:c:v:st:nou:al:oxqkgz:m:ywjUT:OM:K:" opt; do
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
        reaname="$OPTARG"
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
        n_threads=1
        echo "DEPRECATION NOTICE: Option '-k' is deprecated. To disable multi-threading use '-K 1' instead."
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
    U)
        update_manually=true
        ;;
    T)
        user_temp=true
        user_temp_folder=$OPTARG
        ;;
    O)
        force_offline=true
        ;;
    M)
        input_mode=$OPTARG
        ;;
    K)
        n_threads=$OPTARG
        if [ $n_threads -eq 1 ]; then
            use_parallel=false
        fi
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
    #blast_format="qseqid pident evalue bitscore qcovhsp stitle sstart send"
fi

# tmp working directory
fasta=$(readlink -f "$1") # save input file before changing to temporary directory
tmp_fasta=$(basename "${fasta}" .gz | tr ' ' '_')
if [[ "$user_temp" = true ]]; then
    mkdir -p $user_temp_folder
    tmpdir=$(mktemp -d $user_temp_folder/"$tmp_fasta"_XXXXXX)
else
    tmpdir=$(mktemp -d)
    trap 'rm -rf "$tmpdir"' EXIT
fi
echo $tmpdir
cd $tmpdir

# get fasta file
if [[ "$fasta" == *.gz ]]; then # in case fasta is in a archive
    gunzip -c "$fasta" > "$tmp_fasta"
    fasta="$tmp_fasta"
fi
[[ ! -s "$fasta" ]] && { echo Invalid file: $1; exit 0; }
tmpvar=$(basename "$fasta")
fastaID="${tmpvar%.*}"

# Determine if fasta is nucl or prot
if [ $input_mode == "auto" ]; then
    n_char=`cat $fasta | grep -v "^>" | awk '{for(i=1;i<=NF;i++)if(!a[$i]++)print $i}' FS="" | wc -l`
    if [ $n_char -ge "15" ]; then
        echo "Protein fasta detected."
        input_mode="prot"
    else
        echo "Nucleotide fasta detected."
        input_mode="nucl"
    fi
fi

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
    pred_biom=$($dir/predict_biomass_from16S.sh "$fasta")
    if [ "$pred_biom" == "Gram_neg" ] || [ "$pred_biom" == "Gram_pos" ]; then
        taxonomy=Bacteria
    elif [ "$pred_biom" == "Archaea" ]; then
        taxonomy=Archaea
    else
        echo "Taxonomy could be predicted automatically. Assuming default case: Bacteria (Use '-t' parameter to modify it)."
        taxonomy=Bacteria
    fi
    echo Predicted taxonomy: $taxonomy
fi
[[ "$taxonomy" == "bacteria" ]] && taxonomy=Bacteria
[[ "$taxonomy" == "archaea" ]] &&  taxonomy=Archaea

# squence directory
export LC_NUMERIC="en_US.UTF-8"
seqpath=$dir/../dat/seq/$taxonomy
seqpath_user=$dir/../dat/seq/$taxonomy/user
mkdir -p $seqpath/rev $seqpath/unrev $seqpath_user

#check for updates if internet connection is available
if [[ "$force_offline" = false ]]; then
    is_online=$(wget -q --spider http://rz.uni-kiel.de)
    is_running=$(pidof -x "$script_name" -o $$) # do not check for updates if running in parallel mode
    if [[ $is_online -eq 0 && -z "$is_running" ]]; then
        $dir/update_sequences.sh $taxonomy
    fi
    if [[ ! -f $seqpath/rev/sequences.tar.gz  ]] || [[ ! -f $seqpath/unrev/sequences.tar.gz ]] || [[ ! -f $seqpath/rxn/sequences.tar.gz ]]; then
        echo ATTENTION: gapseq sequence archives are missing! Sequences will be needed to be downloaded from uniprot directly which is rather slow.
    fi
fi
download_log=$(mktemp) # remember downloaded files
function already_downloaded(){ 
    if grep -q $1 $download_log; then
        return 0
    else
        return 1
    fi
}


if [ -n "$ecnumber" ] || [ -n "$reaname" ]; then
    # create dummpy pwy template for given ec number
    if [[ -z "$ecnumber" ]]; then
        rea_count=$(echo "$reaname" | tr ';' '\n' | wc -l)
        ecnumber=$(echo "$reaname" | grep -o ";" | tr ';' ',') # get dummy empty comma seperated ec numbers
        pwyname="$reaname"
    elif [[ -z "$reaname" ]]; then
        rea_count=$(echo $ecnumber | tr ',' '\n' | wc -l)
        reaname=$(echo $ecnumber | grep -o "," | tr -d '\n' | tr ',' ';') # get dummy empty colon seperated reaction names
        pwyname=$ecnumber
    else  
        rea_count=$(echo $ecnumber | tr ',' '\n' | wc -l)
    fi
    rea_id=$(seq 1 $rea_count | awk '{print "reaction"$1}' | tr '\n' ',' | sed 's/,$//g')
    pwyKey="custom"
    #pwyDB=$(echo -e "custom\t$ecnumber\t\t\t\t${rea_id::-1}\t$ecnumber\t\t$reaname") # slicing is incompatible?
    pwyDB=$(echo -e "custom\t$pwyname\t\t\t\t${rea_id}\t$ecnumber\t\t"$reaname"")
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
    kegg=$(grep -wFe "$rea" $metaRea | awk -F "\t" {'print $5'})
    altec=""

    for i in "${!ec[@]}"; do
        # 1) search in reaction db by EC
        if [[ -n "${EC_test[i]}" ]]; then
            if [ "$database" == "vmh" ]; then
                dbhit=$(grep -wF ${ec[i]} $reaDB1 | awk -F ',' '{print $1}')
            elif [ "$database" == "seed" ]; then
                dbhit=$(cat $seedEC | cut -f1,3 | grep -wF ${ec[i]} | cut -f1 | tr '|' ' ')
                #dbhit="$dbhit $(grep -wF ${ec[i]} $reaDB4 | awk -F '\t' '{print $4}' | tr '\n' ' ')" # mnxref considers also obsolete seed hits 
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
        if [[ -n "${EC_test[i]}" ]]; then
            altec=$(grep -wF ${ec[i]} $altecdb | grep -P "([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)" -o | grep -vw ${ec[i]})
            altec_src=altec.csv
            if [ -z "$altec" ]; then  
                brendaec=$(grep -wF ${ec[i]} $brenda | grep -P "([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)" -o | grep -vw ${ec[i]})
                [[ `echo "$brendaec" | wc -l` -le 3 ]] && altec=$brendaec # only take unique hits (too many multiple transferred ECs)
                altec_src=brenda
            fi

            if [ "$database" == "vmh" ]; then
                [ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo ${altec//./\\.} | tr ' ' '|')" $reaDB1 | awk -F ',' '{print $1}')" # take care of multiple EC numbers
            elif [ "$database" == "seed" ]; then
                [ -n "$altec" ] && dbhit="$dbhit $(cat $seedEC | cut -f1,3 | grep -wE "$(echo ${altec//./\\.} | tr ' ' '|')" | cut -f1 | tr '|' ' ')" # take care of multiple EC numbers
                #[ -n "$altec" ] && dbhit="$dbhit $(grep -wE "$(echo ${altec//./\\.} | tr ' ' '|')" $reaDB4 | awk -F '\t' '{print $4}')" # mnxref considers also obsolete seed hits 
            fi
        fi
    done
    for aec in $altec; do
        altec_unspecific=$(echo $aec | grep -P "(\\.99\\.[0-9]+$)") # do not accept EC number with unknown acceptors as alternatives
        #echo $aec $altec_unspecific
        if [[ ! " ${ec[@]} " =~ " ${aec} " ]] && [[ -z "$altec_unspecific" ]]; then
            ec+=($aec)
            EC_test+=($(if [[ $aec =~ $re ]]; then echo ${BASH_REMATCH[1]}; fi))
            EC_src+=($altec_src)
        fi
    done

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
if [ "$input_mode" == "nucl" ]; then
    makeblastdb -in "$fasta" -dbtype nucl -out orgdb >/dev/null
fi
if [ "$input_mode" == "prot" ]; then
    makeblastdb -in "$fasta" -dbtype prot -out orgdb >/dev/null
    #diamond makedb -p 16 --in "$fasta" --quiet -d orgdb >/dev/null
fi


cand=""     #list of candidate reactions to be added
bestPwy=""  # list of found pathways
echo -e "ID\tName\tPrediction\tCompleteness\tVagueReactions\tKeyReactions\tKeyReactionsFound\tReactionsFound" > output.tbl # pahtway statistics file

#taxRange=Proteobacteria
if [ -n "$taxRange" ]; then
    validTax=$(grep -i $taxRange $dir/../dat/taxonomy.tbl | cut -f1 | tr '\n' '|' | sed 's/.$//')
    [[ -z "$validTax" ]] && { echo "Taxonomic range not found: $taxRange (available ranges: $dir/../dat/taxonomy.tbl)"; exit 0; }
    pwyDB_new=$(echo "$pwyDB" | grep -wE `echo "TAX-($validTax)"`)
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
        ec_all=$(echo $ecs | awk -v j=$j -F ',' '{print $j}')
        IFS='/' read -r -a ec <<< "$ec_all" # put splitted string into array
        re="([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)"
        EC_test=()
        EC_src=()
        EC_test_bool=false # at least one full EC number found for reaction
        for i in "${!ec[@]}"
        do
            EC_test[i]=$(if [[ ${ec[i]} =~ $re ]]; then echo ${BASH_REMATCH[1]}; fi) # check if not trunked ec number (=> too many hits)
            EC_src[i]=metacyc
            #echo EC check: $i ${ec[i]} ${EC_test[i]}
            [[ -n "{$EC_test[i]}" ]] && EC_test_bool=true
        done
        rea=$(echo $reaids | awk -v j=$j -F ',' '{print $j}')
        reaName=$(echo $reaNames | awk -v j=$j -F ';' '{print $j}' | tr -d '|')
        geneName=$(grep -wFe $rea $metaGenes | awk -F , {'print $2'})
        geneRef=$(grep -wFe $rea $metaGenes | awk -F , {'print $5'})
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
        query=$(mktemp)
        query_all=$(mktemp)
        for i in "${!ec[@]}"; do
            #echo test: ${ec[i]} ${EC_test[i]}
            if [[ -n "${EC_test[i]}" ]]; then
                # check if sequence is not available => try to download
                if [[ ! -f $seqpath/rev/${ec[i]}.fasta || "$update_manually" = true ]]; then
                    if ! already_downloaded "$seqpath/rev/${ec[i]}.fasta"; then
                        [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading reviewed sequences for: ${ec[i]}
                        $dir/uniprot.sh -e "${ec[i]}" -t "$taxonomy" -i $uniprotIdentity -o >/dev/null
                        echo $seqpath/rev/${ec[i]}.fasta >> $download_log
                    fi
                fi
                if [[ (! -f $seqpath/unrev/${ec[i]}.fasta && $seqSrc -gt 1) || "$update_manually" = true ]]; then 
                    if ! already_downloaded "$seqpath/unrev/${ec[i]}.fasta"; then
                        [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading unreviewed sequences for: ${ec[i]}
                         $dir/uniprot.sh -u -e "${ec[i]}" -t "$taxonomy" -i $uniprotIdentity -o >/dev/null
                        echo $seqpath/unrev/${ec[i]}.fasta >> $download_log
                    fi
                fi
                if [ -s "$seqpath_user/${ec[i]}.fasta" ]; then
                    query_tmp=$seqpath_user/${ec[i]}.fasta
                elif [ $seqSrc -eq 1 ]; then
                    query_tmp=$seqpath/rev/${ec[i]}.fasta
                elif [ $seqSrc -eq 2 ]; then
                    query_tmp=$seqpath/rev/${ec[i]}.fasta
                    [[ ! -s $query_tmp ]] && query_tmp=$seqpath/unrev/${ec[i]}.fasta
                elif [ $seqSrc -eq 3 ]; then
                    query_tmp=$(mktemp)
                    cat $seqpath/rev/${ec[i]}.fasta $seqpath/unrev/${ec[i]}.fasta | awk '/^>/{f=!d[$1];d[$1]=1}f' > $query_tmp # use awk to remove duplicates
                elif [ $seqSrc -eq 4 ]; then
                    query_tmp=$seqpath/unrev/${ec[i]}.fasta
                fi
                if [[ -s $query_tmp ]]; then
                    cat $query_tmp >> $query_all
                    if [[ "$query_tmp" == "$seqpath_user"* ]]; then
                        [[ verbose -ge 1 ]] && echo -e "\t\t--> Found user sequences: $query_tmp (`cat $query_tmp | grep ">" | wc -l` sequences)" 
                    elif [[ $i -eq 0 ]]; then
                        [[ verbose -ge 1 ]] && echo -e "\t\t--> Found sequences: $query_tmp (`cat $query_tmp | grep ">" | wc -l` sequences)"
                    else
                        [[ verbose -ge 1 ]] && echo -e "\t\t--> Found alternative EC ${ec[i]} from ${EC_src[i]}: $query_tmp (`cat $query_tmp | grep ">" | wc -l` sequences)"
                    fi
                fi
            fi
        done
        [[ -s $query_all ]] && { cat $query_all | awk '/^>/{f=!d[$1];d[$1]=1}f' > $query; } # no duplicates
        ec_avail=$(join_by / "${EC_test[@]}") # all valid ec numbers
        
        # if no EC number is available or no sequence was found for EC number then use reaction name instead for sequence search
        if [[ -n "$reaName" ]] && ( [[ "$EC_test_bool" = false ]] || [[ ! -s "$query" ]] );then
            reaNameHash=$(echo -n "$reaName" | md5sum | awk '{print $1}')
            # check if sequence is not available => try to download
            if [[ ! -f $seqpath/rev/$reaNameHash.fasta  || "$update_manually" = true ]]; then
                if ! already_downloaded "$seqpath/rev/$reaNameHash.fasta"; then
                    [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading reviewed sequences for: $reaName "\n\t\t(hash: $reaNameHash)" 
                    $dir/uniprot.sh -r "$reaName" -t "$taxonomy" -i $uniprotIdentity -o >/dev/null
                    echo $seqpath/rev/$reaNameHash.fasta >> $download_log
                fi
            fi
            if [[ (! -f $seqpath/unrev/$reaNameHash.fasta && $seqSrc -gt 1) || "$update_manually" = true ]]; then 
                if ! already_downloaded "$seqpath/unrev/$reaNameHash.fasta"; then
                    [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading unreviewed sequences for: $reaName "\n\t\t(hash: $reaNameHash)" 
                    $dir/uniprot.sh -u -r "$reaName" -t "$taxonomy" -i $uniprotIdentity -o >/dev/null
                    echo $seqpath/unrev/$reaNameHash.fasta >> $download_log
                fi
            fi
            if [ -s "$seqpath_user/$reaNameHash.fasta" ]; then
                query=$seqpath_user/$reaNameHash.fasta
            elif [ $seqSrc -eq 1 ]; then
                query=$seqpath/rev/$reaNameHash.fasta
            elif [ $seqSrc -eq 2 ]; then
                query=$seqpath/rev/$reaNameHash.fasta
                [[ ! -s $query ]] && query=$seqpath/unrev/$reaNameHash.fasta
            elif [ $seqSrc -eq 3 ]; then
                query=$(mktemp)
                cat $seqpath/rev/$reaNameHash.fasta $seqpath/unrev/$reaNameHash.fasta | awk '/^>/{f=!d[$1];d[$1]=1}f' > $query
            elif [ $seqSrc -eq 4 ]; then
                query=$seqpath/unrev/$reaNameHash.fasta
            fi
            if [[ -s $query ]]; then
                if [[ "$query" == "$seqpath_user"* ]]; then
                    [[ verbose -ge 1 ]] && echo -e "\t\t--> Found user sequences: $query (`cat $query | grep ">" | wc -l` sequences)"
                else
                    [[ verbose -ge 1 ]] && echo -e "\t\t--> Found sequences: $query (`cat $query | grep ">" | wc -l` sequences)"
                fi
            fi
        fi

        # sequence by gene name
        if [[ -n "$geneName" ]] && [[ -n "$geneRef" ]] && [[ "$use_gene_seq" = true ]]; then
            if [[ ! -f $seqpath/rxn/$rea.fasta || "$update_manually" = true ]]; then
                reaSeqTmp=$(mktemp)
                for gr in $geneRef
                do
                    if ! already_downloaded "$seqpath/rxn/$gr.fasta"; then
                        [[ verbose -ge 1 ]] && echo -e '\t\t'Downloading sequences for: $gr
                        $dir/uniprot.sh -d $gr -t "$taxonomy" -i $uniprotIdentity -o >/dev/null
                        echo $seqpath/rxn/$gr.fasta >> $download_log
                    fi
                    if [ -f $seqpath/rxn/$gr.fasta ]; then
                        cat $seqpath/rxn/$gr.fasta >> $reaSeqTmp
                        # rm $seqpath/rxn/$gr.fasta # shouldn't be deleted to allow monitoring of changes 
                    fi
                done
                if [ -s "$reaSeqTmp" ]; then
                   mv $reaSeqTmp $seqpath/rxn/$rea.fasta
                else
                    touch $seqpath/rxn/$rea.fasta # create empty file if no gene seq data is found to avoid reoccuring download attempt
                fi
            fi
            
            if [ -s "$seqpath_user/$rea.fasta" ]; then
                query_gene=$seqpath_user/$rea.fasta
            else
                query_gene=$seqpath/rxn/$rea.fasta
            fi
            #merge sequence data
            if [[ -s $query_gene ]]; then
                if [[ "$query_gene" == "$seqpath_user"* ]]; then
                    [[ verbose -ge 1 ]] && echo -e "\t\t--> Found user sequences: $query_gene (`cat $query_gene | grep ">" | wc -l` sequences)"
                else
                    [[ verbose -ge 1 ]] && echo -e "\t\t--> Found sequences: $query_gene (`cat $query_gene | grep ">" | wc -l` sequences)"
                fi
                query_merge2=$(mktemp)
                cat $query $query_gene | awk '/^>/{f=!d[$1];d[$1]=1}f' > $query_merge2 # no duplicates
                query=$query_merge2
            fi
        
        fi
        
        [[ "$skipBlast" = true ]] && { continue; }

        if [ -s $query ]; then
            [[ verbose -ge 1 ]] && echo -e "\t\tFinal file: $query (`cat $query | grep ">" | wc -l` sequences)"
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
                        if { { ! [ -x "$(command -v parallel)" ]; } && [ "$input_mode" == "nucl" ]; } || [ "$use_parallel" = false ]; then # try to use parallelized version # TODO CORRECT THIS!
                            if [ "$input_mode" == "nucl" ]; then
                                #echo "Hello. We are nucl here! (ncpu=1)"
                                tblastn -db orgdb -query $q -qcov_hsp_perc $covcutoff -outfmt "6 $blast_format" > query.blast
                            fi
                            if [ "$input_mode" == "prot" ]; then
                                #echo "Hello. We are proteins here!"
                                blastp -db orgdb -query $q -qcov_hsp_perc $covcutoff -outfmt "6 $blast_format" > query.blast
                            fi
                        else
                            if [ "$input_mode" == "nucl" ]; then
                                cat $q | parallel --gnu --will-cite --block 50k -j $n_threads --recstart '>' --pipe tblastn -db orgdb -qcov_hsp_perc $covcutoff -outfmt \'"6 $blast_format"\' -query - > query.blast
                                #tblastn -db orgdb -query $q -qcov_hsp_perc $covcutoff -num_threads $n_threads -outfmt "6 $blast_format" > query.blast
                            fi
                            if [ "$input_mode" == "prot" ]; then
                                #cat $q | parallel --gnu --will-cite --block 50k --recstart '>' --pipe blastp -db orgdb -qcov_hsp_perc $covcutoff -num_threads 16 -outfmt \'"6 $blast_format"\' -query - > query.blast
                                blastp -db orgdb -query $q -qcov_hsp_perc $covcutoff -num_threads $n_threads -outfmt "6 $blast_format" > query.blast
                                #diamond blastp -d orgdb -q $q -b6 --query-cover $covcutoff --outfmt 6 $blast_format -p $n_threads > query.blast
                            fi
                        fi
                        cat query.blast >> $out
                        bhit=$(cat query.blast | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff_tmp -v covcutoff=$covcutoff '{if ($2>=identcutoff && $4>=bitcutoff && $5>=covcutoff) print $0}')
                        if [ -n "$bhit" ]; then
                            bestsubunithit=$(echo "$bhit" | sort -rgk 4,4)
                            [[ verbose -ge 1 ]] && [[ $iterations -gt 1 ]] && echo "$bestsubunithit" | head -1 | cut -f1 | grep -f - $query | sed "s/^/\t\t\t$subunit_id hit: /" 
                            ((subunits_found++))
                            [[ "$subunit_id" == "Subunit undefined" ]] && ((subunits_undefined_found++))
                            [[ $iterations -gt 1 ]] && echo "$bestsubunithit" | awk -v exception="$is_exception" -v subunit="$subunit_id" -v rea="$rea" -v reaName="$reaName" -v ec=$ec_avail -v dbhit="$dbhit" -v pwy="$pwy" '{print rea"\t"reaName"\t"ec"\t"NA"\t"$0"\t"pwy"\t""good_blast""\t""NA""\t"dbhit"\t"subunit"\t"exception"\t""NA"}' >> reactions.tbl
                            [[ "$exhaustive" = false ]] && break 
                        fi
                        rm query.blast
                    done
                    if [[ $subunits_found -eq $subunits_found_old ]] && [[ $iterations -gt 1 ]]; then
                        if [ "$includeSeq" = true ]; then
                            echo -e "$rea\t$reaName\t$ec_avail\tNA\t\t\t\t\t\t\t\t\t\t$pwy\tno_blast\tNA\t$dbhit\t$subunit_id\t$is_exception\tNA" >> reactions.tbl # subunit not found 
                        else
                            echo -e "$rea\t$reaName\t$ec_avail\tNA\t\t\t\t\t\t\t\t\t$pwy\tno_blast\tNA\t$dbhit\t$subunit_id\t$is_exception\tNA" >> reactions.tbl # subunit not found 
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
                    [[ $iterations -le 1 ]] && echo "$besthit_all" | awk -v exception="$is_exception" -v rea="$rea" -v reaName="$reaName" -v ec=$ec_avail -v dbhit="$dbhit" -v pwy="$pwy" '{print rea"\t"reaName"\t"ec"\t"NA"\t"$0"\t"pwy"\t""good_blast""\t""NA""\t"dbhit"\t""NA""\t"exception"\t""NA"}' >> reactions.tbl # only for non-subunit hits
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
                    
                    if [[ -n "$dbhit" && "$dbhit" != " " ]]; then
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
                        echo -e "$rea\t$reaName\t$ec_avail\tNA\t$somehit_all\t$pwy\tbad_blast\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
                    else
                        [[ verbose -ge 1 ]] && echo -e '\t\t'NO hit because of missing subunits
                        if [[ "$subunits_former_run" = true ]];then # log also subunits from former run
                            tmp_log=$(cat reactions.tbl | awk -F '\t' -v rea="$rea" -v reaName="$reaName" -v ec="$ec_avail" -v pwy="$pwy" '{OFS=FS} $1==rea && $2==reaName && $3==ec {$13=pwy; print}')
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
                        echo -e "$rea\t$reaName\t$ec_avail\tNA\t\t\t\t\t\t\t\t\t\t$pwy\tno_blast\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
                    else
                        echo -e "$rea\t$reaName\t$ec_avail\tNA\t\t\t\t\t\t\t\t\t$pwy\tno_blast\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
                    fi
                fi
            fi
        else
            [[ verbose -ge 1 ]] && echo -e "\t\tNo sequence data found"
            ((vague++))
            [[ -n "$dbhit" ]] && pwyVage="$pwyVage$dbhit "
            [[ $keyRea = *"$rea"* ]] && vagueKeyReaFound="$vagueKeyReaFound $rea"
            if [ "$includeSeq" = true ]; then
                echo -e "$rea\t$reaName\t$ec_avail\tNA\t\t\t\t\t\t\t\t\t\t$pwy\tno_seq_data\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
            else
                echo -e "$rea\t$reaName\t$ec_avail\tNA\t\t\t\t\t\t\t\t\t$pwy\tno_seq_data\tNA\t$dbhit\tNA\t$is_exception\tNA" >> reactions.tbl 
            fi
        fi
        [[ verbose -ge 2 ]] && echo -e "\t\tCandidate reactions: $dbhit"
    done # pathway
    
    if [ $count -eq 0 ]; then
        completeness=0
    else
        check_vague=$(echo "$vague < $count*$vagueCutoff" | bc) # vague reactions shouldn't make more than certain amount of total reactions
        if [ "$strictCandidates" = false ] && [ $check_vague -eq 1 ] ; then # if vague reaction are considered they should not influence the completeness threshold
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
    # add reactions of pathways (even if no blast hit) if above threshold (and no key enzyme is missed)
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
    # B) Consider as complete pathway because of completeness threshold (key enzymes should be present too)
    elif [[ $completeness -ge $completenessCutoffNoHints ]] && [[ "$KeyReaFracAvail" -eq 1 ]] && [[ "$strictCandidates" = false ]]; then
        pwy_status="threshold"
        [[ verbose -ge 1 ]] && echo "Consider pathway to be present because of completeness threshold!"
        prediction=true
        cand="$cand$pwyVage$pwyNoHitFound "
        bestPwy="$bestPwy$name ($completeness% completeness, added because of completeness threshhold)\n"
    # C) Consider as complete pathway because of key enzymes (lower threshold)
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

# add gapseq version and sequence database status to table comments head
gapseq_version=`$dir/.././gapseq -v`
seqdb_version=`md5sum $dir/../dat/seq/$taxonomy/rev/sequences.tar.gz | cut -c1-7`
seqdb_date=$(stat -c %y $dir/../dat/seq/$taxonomy/rev/sequences.tar.gz | cut -c1-10)

sed -i "1s/^/# $gapseq_version\n/" $curdir/${fastaID}-$output_suffix-Reactions.tbl
sed -i "2s/^/# Sequence DB md5sum: $seqdb_version ($seqdb_date, $taxonomy)\n/" $curdir/${fastaID}-$output_suffix-Reactions.tbl
sed -i "1s/^/# $gapseq_version\n/" $curdir/${fastaID}-$output_suffix-Pathways.tbl
sed -i "2s/^/# Sequence DB md5sum: $seqdb_version ($seqdb_date, $taxonomy)\n/" $curdir/${fastaID}-$output_suffix-Pathways.tbl

# print annotation genome coverage
[[ verbose -ge 1 ]] && [[ "$anno_genome_cov" = true ]] && Rscript $dir/coverage.R "$fasta" $curdir/${fastaID}-$output_suffix-Reactions.tbl 

# cleaning
[[ -s "$tmp_fasta" ]] && rm "$tmp_fasta"


ps -p $$ -o %cpu,%mem,cmd
end_time=`date +%s`
echo Running time: `expr $end_time - $start_time` s.
