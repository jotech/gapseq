#!/bin/bash

start_time=`date +%s`

pathways=""
database="seed"
pwyDatabase="metacyc,custom"
verbose=1
taxonomy="Bacteria"
taxRange="all" # taxonomic range for pathways
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
exhaustive=false
seqSrc=2
anno_genome_cov=false
use_gene_seq=true
stop_on_files_exist=false
update_manually=false
user_temp=false
gramstaining=NA
force_offline=false
input_mode="auto"
output_dir=.
aliTool="blast"
OS=$(uname -s)
if [ "$OS" = "Darwin" -o "$OS" = "FreeBSD" ]; then
    n_threads=$(sysctl hw.ncpu|cut -f2 -d' ')
else
    n_threads=`grep -c ^processor /proc/cpuinfo`
fi

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
    echo "  -f Path to directory, where output files will be saved (default: current directory)"
    echo "  -U Do not use gapseq sequence archive and update sequences from uniprot manually (very slow) (default: $update_manually)"
    echo "  -T Set user-defined temporary folder (default: $user_temp)"
    echo "  -O Force offline mode (default: $force_offline)"
    echo "  -M Input genome mode. Either 'nucl', 'prot', or 'auto' (default '$input_mode')"
    echo "  -K Number of threads for sequence alignments. If option is not provided, number of available CPUs will be automatically determined."
    echo "  -A Tool to be used for sequence alignments (blast, mmseqs2, diamond; default: $aliTool)"
    echo ""
    echo "Details:"
    echo "\"-t\": if 'auto', gapseq tries to predict if the organism is Bacteria or Archaea based on the provided genome sequence. The prediction is based on the 16S rRNA gene sequence using a classifier that was trained on 16S rRNA genes from organisms with known Gram-staining phenotype. In case no 16S rRNA gene was found, a k-mer based classifier is used instead."

exit 1
}


# paths and variables
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
script_name=$(basename -- "$0")
uniprotIdentity=0.9 # clustered uniprot database (0.5 or 0.9) # This variable has currently not effect. Reviewed cluster identity: 0.9, unreviews: 0.5. These values are hardcoded in 'src/uniprot.sh'
metaPwy=$dir/../dat/meta_pwy.tbl
keggPwy=$dir/../dat/kegg_pwy.tbl
seedPwy=$dir/../dat/seed_pwy.tbl
customPwy=$dir/../dat/custom_pwy.tbl
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
metaGenes=$dir/../dat/meta_genes.csv

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "h?p:e:r:d:i:b:c:v:st:nou:al:oxqkgz:m:ywjf:UT:OM:K:A:" opt; do
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
    f)
        output_dir=$OPTARG
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
        ;;
    A)
        aliTool=$(echo "$OPTARG" | tr '[:upper:]' '[:lower:]')
        if [[ "$aliTool" != "blast" && "$aliTool" != "diamond" && "$aliTool" != "mmseqs2" ]]; then
            echo "Error: Invalid value for -A. Expected 'blast', 'diamond', or 'mmseqs2', got '$OPTARG'."
            exit 1
        fi
        ;;
    esac
done
shift $((OPTIND-1))
[ "$1" = "--" ] && shift

# after parsing arguments, only fasta file should be there
[ "$#" -ne 1 ] && { usage; }

# alignment statistics format
if [ "$includeSeq" = true ]; then
    blast_format="qseqid pident evalue bitscore qcovs stitle sstart send sseq"
    mmseqs_format="qheader,pident,evalue,bits,qcov,theader,tstart,tend,qseq"
    diamnd_format="qseqid pident evalue bitscore qcovhsp stitle sstart send sseq"
else
    blast_format="qseqid pident evalue bitscore qcovs stitle sstart send"
    mmseqs_format="qheader,pident,evalue,bits,qcov,theader,tstart,tend"
    diamond_format="qseqid pident evalue bitscore qcovhsp stitle sstart send"
fi

# set output directory
case $output_dir in
    /*)
        # absolute path
        output_dir=$output_dir
        ;;
    ~*)
        # relative to $HOME directory
        output_dir="${output_dir/#\~/$HOME}"
        ;;
    *)
        # relative path to current directory
        output_dir=$curdir/$output_dir
        ;;
esac

# create path if it does not yet exist
mkdir -p $output_dir || { echo "$output_dir not writable. Aborting..."; exit 1; }

# tmp working directory
fasta=$(readlink -f "$1") # save input file before changing to temporary directory
tmp_fasta=$(basename "${fasta}" .gz | tr ' ' '_')
if [[ "$user_temp" = true ]]; then
    mkdir -p $user_temp_folder || { echo "Temporary directory $user_temp_folder cannot be created. Aborting..."; exit 1; }
    tmpdir=$(mktemp -d $user_temp_folder/"$tmp_fasta"_XXXXXX)
    [[ "$tmpdir" != /* ]] && tmpdir="$(cd "$tmpdir" && pwd)" # if relative path — convert it to absolute
else
    tmpdir=$(mktemp -d)
fi
trap 'rm -rf "$tmpdir"' EXIT
[[ $verbose -ge 1 ]] && echo $tmpdir
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
    input_mode=`$dir/./nuclprot.sh $fasta`
    
    if [ $input_mode == "prot" ]; then
        [[ $verbose -ge 1 ]] && echo "Protein fasta detected."
    else
        [[ $verbose -ge 1 ]] && echo "Nucleotide fasta detected."
    fi
    
fi

if [ $input_mode == "nucl" ]; then
    newtranslate=true
    # Check if genome was already translated
    if [ -f $output_dir/${fastaID}.faa.gz ]; then
        # Check if contigs of found ORFs matches contig names in nucleotide fasta
        faacont=`zcat $output_dir/${fastaID}.faa.gz | grep "^>" | sed -E 's/^>(.+)_[0-9]+ # .*/\1/' | sort -u`
        fnacont=`cat $fasta | grep "^>" | sed -E 's/^>([^ ]+).*/\1/' | sort -u`
        reusefaa=true
        for entry in $faacont; do
            if ! echo "$fnacont" | grep -qx "$entry"; then
                reusefaa=false
                break
            fi
        done
        if [ $reusefaa == "true" ]; then
            [[ $verbose -ge 1 ]] && echo "Re-using previously translated genome: $output_dir/${fastaID}.faa.gz"
            gunzip -c "$output_dir/${fastaID}.faa.gz" > "$fastaID.faa"
            fasta="$fastaID.faa"
            newtranslate=false
        fi
    fi
    
    if [ $newtranslate == "true" ]; then
        [[ $verbose -ge 1 ]] && echo -n "Translating genomic nucleotide fasta to protein fasta..."
        $dir/translate_genome.sh -i "$fasta" -o "$fastaID" -K $n_threads 
        fasta="$fastaID.faa"
        transl_table=`cat ${fastaID}_code`
        rm ${fastaID}_code
        orf_count=$(grep -c "^>" "$fasta")
        [[ $verbose -ge 1 ]] && echo "$orf_count ORFs (translation table: $transl_table)"
    fi
fi

# pathways or ec number as well as fasta file have to be provided
( [ -z "$pathways" ] && [ -z "$ecnumber" ] && [ -z "$reaname" ] ) || [ -z "$fasta" ]  && { usage; }

# select pathway keys to be used in database search
case $pathways in
    all)
        pwyKey="Pathways|Enzyme-Test|seed|kegg"
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
    min)
        pwyKey="\\|ETOH-ACETYLCOA-ANA-PWY\\||\\|GLNSYN-PWY\\||\\|GLUCONEO-PWY\\||\\|GLUGLNSYN-PWY\\||\\|GLUTAMATE-DEG1-PWY\\||\\|GLUTAMATE-SYN2-PWY\\||\\|GLUTAMINEFUM-PWY\\||\\|GLUTSYNIII-PWY\\||\\|GLYCOLYSIS\\||\\|GLYOXYLATE-BYPASS\\||\\|NONOXIPENT-PWY\\||\\|OXIDATIVEPENT-PWY\\||\\|P185-PWY\\||\\|P21-PWY\\||\\|PWY0-1312\\||\\|PWY0-1315\\||\\|PWY0-1329\\||\\|PWY0-1334\\||\\|PWY0-1335\\||\\|PWY0-1353\\||\\|PWY0-1517\\||\\|PWY0-1565\\||\\|PWY0-1567\\||\\|PWY0-1568\\||\\|PWY-4341\\||\\|PWY-5084\\||\\|PWY-5480\\||\\|PWY-5482\\||\\|PWY-5484\\||\\|PWY-5690\\||\\|PWY-5766\\||\\|PWY-5913\\||\\|PWY-6028\\||\\|PWY-6333\\||\\|PWY-6543\\||\\|PWY-6549\\||\\|PWY66-21\\||\\|PWY66-398\\||\\|PWY-6697\\||\\|PWY-6964\\||\\|PWY-7167\\||\\|PWY-7685\\||\\|PWY-7686\\||\\|PWY-7980\\||\\|PWY-8178\\||\\|PWY-8215\\||\\|PWY-8274\\||\\|PWY-8404\\||\\|PYRUVDEHYD-PWY\\||\\|TCA-1\\||\\|TCA\\|"
        ;;
    small)
        tca="\\|TCA\\||P105-PWY"
        resp="PWY-7544|PWY0-1334"
        ferm="PYR-TO-BUT-NADPH|2FLHYD|PWY-7385|PWY-6344|PWY-5938|PWY-5494|PWY-5497|P164-PWY|PWY-8086|PWY-5677|PWY-8014|P162-PWY|CENTFERM-PWY|P122-PWY|PWY-6130|P108-PWY|PYR-ACCOA-FROX|BIFIDOSHUNT2|LNT-DEGRADATION|PWY-6583|PWY-5437|FERMENTATION-PWY"
        glyc="ANAGLYCOLYSIS-PWY|PWY-1042"
        ppp="NONOXIPENT-PWY|OXIDATIVEPENT-PWY"
        aasyn="PWY-4361|PWY-3941|ARGASEDEG-PWY|ASPARAGINESYN-PWY|PWY-2942|ALANINE-VALINESYN-PWY|PHESYN|PWY-2941|TRPSYN-PWY2|PWY-4341|SERSYN-PWY|GLYSYN-THR-PWY|ASPARAGINE-BIOSYNTHESIS|GLUTSYNIII-PWY|TRPSYN-PWY|PWY-3461|PWY-5101|LEUSYN-PWY|PWY-3081|PWY-5097|ILEUSYN-PWY|GLUTAMATE-SYN2-PWY|HOMOSERSYN-PWY|PWY-3982|PROSYN-PWY|PWY-5344|HOMOSER-THRESYN-PWY|GLNSYN-PWY|HOMOCYSDEGR-PWY|CYSTSYN-PWY|HISTSYN-PWY|GLYSYN-ALA-PWY|HOMOSER-METSYN-PWY|ARGSYNBSUB-PWY|PWY0-901|PWY490-4|DAPLYSINESYN-PWY"
        cosyn="PWY-6148|PWY-8171|PWY0-1275|PWY-5443|MENAQUINONESYN-PWY|PWY-6892|PWY-6383|PWY-7761|PWY-7165|NADPHOS-DEPHOS-PWY|PWY-7377|PWY0-1507|PWY-5509|PWY-5837|PWY-6893|PWY-7356|NAD-BIOSYNTHESIS-II|PWY-7376|PWY-6151|RIBOSYN2-PWY|PWY-7962|PWY-7729|NADPHOS-DEPHOS-PWY-1|PWY-2161|PWY-6910|PYRIDNUCSYN-PWY|NADP-RED-FD-NADH|PWY-6890|PANTO-PWY|PWY-7378|P241-PWY|COA-PWY|PWY-6614|PWY-5653|PYRIDNUCSAL-PWY|PWY-5785|PYRIDOXSYN-PWY|PWY-6891|1CMET2-PWY"
        nusyn="PWY-7199|IMPSYN-HIGH-HCO3|PWY-7193|PWY-6610|SALVPURINE2-PWY|PWY-7205|PWY-6599|PWY-7791|PWY-7790|PWY-7187|PWY-7210|P121-PWY|PWY-6545|PWY-7183|PWY-6121|PWY-5686|PWY-6556|PWY-6123|PWY-7221|PWY-7219|PWY-7197|PWY-7176"
        carbodeg="PWY-6572|RIBOKIN-PWY|PWY-7077|PWY-7581|PWY-7247|PWY-5941|2FLHYD|PWY-6717|PWY-6507|PWY0-1314|ARABCAT-PWY|PWY-7242|RHAMCAT-PWY|LACTOSECAT-PWY|DARABITOLUTIL-PWY|PWY0-1309|TREDEGLOW-PWY|PWY0-1324|MANNCAT-PWY|PWY0-1301|PWY-7178|GLYCEROLMETAB-PWY|PWY-2221|PWY-7180|GALACTITOLCAT-PWY|PWY-621|GALACTUROCAT-PWY|PWY-3801|LACTOSEUTIL-PWY|PWY-5517|PWY-8121|PWY-4261|PWY-6130|DARABCATK12-PWY|GLUAMCAT-PWY|BIFIDOSHUNT2|EXC-STARCH-N27|LNT-DEGRADATION|GLUCOSE1PMETAB-PWY|GLYCOCAT-PWY|PWY-6906|PWY-6317"
        carboxdeg="GLUCONSUPER-PWY|IDNCAT-PWY|PWY-7247|PWY-7948|2FLHYD|PWY0-1313|PROPIONMET-PWY|PWY0-42|PWY-7242|PWY-8134|PWY-5177|GALACTARDEG-PWY|GALACTUROCAT-PWY|GLYCOLATEMET-PWY|PWY-6518|GALACTCAT-PWY|BIFIDOSHUNT2|PWY-7754|LNT-DEGRADATION|PWY-6697"
        pwyKey="$tca|$resp|$ferm|$glyc|$ppp|$aasyn|$cosyn|$nusyn|$carbodeg|$carboxdeg"
        noSuperpathways=false # some listed pathways are superpathways
        ;;
    kegg)
        pwyKey=kegg
        ;;
    *)
        pwyKey=$pathways
        ;;
esac

# predict taxonomy (Bacteria or Archaea?)
if [ $taxonomy == "auto" ]; then
    cp $dir/../dat/seq/hmm/domain.hmm.gz .
    gunzip domain.hmm.gz
    hmmsearch --tblout $fastaID.tblout --cpu $n_threads domain.hmm $fasta > /dev/null
    taxonomy=`Rscript $dir/predict_domain.R "$dir" "$fastaID.tblout"`
    rm domain.hmm
    rm $fastaID.tblout
    
    [[ $verbose -ge 1 ]] && echo Predicted taxonomy: $taxonomy
    

fi
[[ "$taxonomy" == "bacteria" ]] && taxonomy=Bacteria
[[ "$taxonomy" == "archaea" ]] &&  taxonomy=Archaea


# Predict Gram staining
if [ "$taxonomy" == "Bacteria" ]; then
    cp $dir/../dat/seq/hmm/gram.hmm.gz .
    gunzip gram.hmm.gz
    hmmsearch --tblout $fastaID.tblout --cpu $n_threads gram.hmm $fasta > /dev/null
    gramstaining=`Rscript $dir/predict_gramstaining.R "$dir" "$fastaID.tblout"`
    rm gram.hmm
    rm $fastaID.tblout

    [[ $verbose -ge 1 ]] && echo Predicted Gram-staining: $gramstaining

fi



# Follow taxonomy prediction for pathway tax range if set to "auto"
if [ "$taxRange" == "auto" ]; then
    taxRange=$taxonomy
fi


# squence directory
export LC_NUMERIC="en_US.UTF-8"
seqpath=$dir/../dat/seq/$taxonomy
seqpath_user=$dir/../dat/seq/$taxonomy/user
mkdir -p $seqpath/rev $seqpath/unrev $seqpath_user

#check for updates if internet connection is available
if [[ "$force_offline" = false ]]; then
    wget -q --spider https://zenodo.org
    is_online=$?
    [[ `pgrep -f $0` != "$$" ]] && is_running=yes
    if [[ $is_online -eq 0 && -z "$is_running" ]]; then
        $dir/update_sequences.sh $taxonomy
    fi
    if [[ ! -f $seqpath/rev/sequences.tar.gz  ]] || [[ ! -f $seqpath/unrev/sequences.tar.gz ]] || [[ ! -f $seqpath/rxn/sequences.tar.gz ]]; then
        echo ATTENTION: gapseq sequence archives are missing! Sequences will be needed to be downloaded from uniprot directly which is rather slow.
    fi
fi
download_log=$(mktemp -p $tmpdir) # remember downloaded files
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
        metaRea_hit=$(grep -wFe "$reaname" $metaRea)
        if [[ -n "$metaRea_hit" ]]; then
            rea_count=1
            ecnumber=$(echo "$metaRea_hit" | awk -F "\t" {'print $4'} | sed 's/EC-//g' | sed 's/,/\//g')
            pwyname="$reaname"
        else
            rea_count=$(echo "$reaname" | tr ';' '\n' | wc -l)
            ecnumber=$(echo "$reaname" | grep -o ";" | tr ';' ',') # get dummy empty comma seperated ec numbers
            pwyname="$reaname"
        fi
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
    [[ $verbose -ge 1 ]] && { echo $pwyDatabase; }
    [[ "$pwyDatabase" =~ "all" ]]     && cat $metaPwy $keggPwy $seedPwy $customPwy > allPwy
    [[ "$pwyDatabase" =~ "metacyc" ]] && cat $metaPwy >> allPwy
    [[ "$pwyDatabase" =~ "kegg" ]]    && cat $keggPwy >> allPwy
    [[ "$pwyDatabase" =~ "seed" ]]    && cat $seedPwy >> allPwy
    [[ "$pwyDatabase" =~ "custom" ]]  && cat $customPwy >> allPwy
    dupli=$(cat allPwy | cut -f1 | sort | uniq -d | tr -d 'id' | sed '/^$/d')
    if [ -n "$dupli" ]; then
        [[ $verbose -ge 1 ]] && echo Duplicated pathway IDs found: $dupli will only use $customPwy
        dupli_search=$(echo "$dupli" | sed 's/|/\\|/g' |tr '\n' '|' | rev | cut -c2- | rev)
        [[ $verbose -ge 1 ]] && echo "$dupli_search"
        cat allPwy | grep -wEv "$dupli_search" > allPwy.tmp
        cat $customPwy | grep -wE "$dupli_search" >> allPwy.tmp
        mv allPwy.tmp allPwy
        #cat allPwy | grep -wE "$dupli_search"
    fi
    # cat allPwy | grep -wEi $pwyKey | wc -l
    pwyDB=$(cat allPwy | grep -wEi $pwyKey | awk -F "\t" '{if ($6) print $0;}')
    NApwy=$(cat allPwy | grep -wEi $pwyKey | awk -F "\t" '{if (!$6) print $1, $2;}')
    [[ -n "$NApwy" ]] && echo Pathways ignored because no reactions found: $(echo "$NApwy" | wc -l)
    [[ -n "$NApwy" ]] && [[ $verbose -ge 2 ]] && echo "$NApwy"
    [[ "$noSuperpathways" = true ]] && pwyDB=$(echo "$pwyDB" | grep -v 'Super-Pathways')
    [ -z "$ecnumber" ] && [ -z "$pwyDB" ] && { echo "No pathways found for key $pwyKey"; exit 1; }
fi
[ -z "$output_suffix" ] && output_suffix=$pathways


# check to quit if output files already exist
if [[ "$stop_on_files_exist" = true ]] && [[ -f $output_dir/${fastaID}-$output_suffix-Pathways.tbl ]] && [[ -f $output_dir/${fastaID}-$output_suffix-Reactions.tbl ]]; then
    echo Pathway and reaction output files already exist, exiting...
    exit 0
fi


# function to get database hits for ec number
getDBhit(){
    kegg=$(grep -wFe "$rea" $metaRea | awk -F "\t" {'print $5'})
    altec=""
    dbhit=""

    for i in "${!ec[@]}"; do
        # 1) search in reaction db by EC
        if [[ -n "${EC_test[i]}" ]]; then
            if [ "$database" == "vmh" ]; then
                dbhit="$dbhit $(grep -wF ${ec[i]} $reaDB1 | awk -F ',' '{print $1}')"
            elif [ "$database" == "seed" ]; then
                dbhit="$dbhit $(cat $seedEC | cut -f1,3 | grep -wF ${ec[i]} | cut -f1 | tr '|' ' ')"
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

    dbhit=$(echo "$dbhit" | sed 's/^[[:space:]]*//')
}



cand=""     #list of candidate reactions to be added
bestPwy=""  # list of found pathways
echo -e "ID\tName\tPrediction\tCompleteness\tVagueReactions\tKeyReactions\tKeyReactionsFound\tReactionsFound" > output.tbl # pahtway statistics file

#taxRange=Proteobacteria
if [ -n "$taxRange" ] && [ "$taxRange" != "all" ] && [ "$pathways" != "custom" ]; then
    validTax=$(grep -i $taxRange $dir/../dat/taxonomy.tbl | cut -f1 | tr '\n' '|' | sed 's/.$//')
    [[ -z "$validTax" ]] && { echo "Taxonomic range not found: $taxRange (available ranges: $dir/../dat/taxonomy.tbl)"; exit 0; }
    pwyDB_new=$(echo "$pwyDB" | grep -wE `echo "TAX-($validTax)"`)
    pwyDB_old=$(echo "$pwyDB" | awk -F '\t' 'BEGIN {OFS=FS="\t"} $5=="" {print $0}')
    pwyDB=$(echo -e "$pwyDB_new\n$pwyDB_old")
    [[ -z "$pwyDB" ]] && { echo "No pathways found"; exit 0; }
fi

pwyNr=$(echo "$pwyDB" | wc -l)
[[ $verbose -ge 1 ]] && echo Checking for pathways and reactions in: $1 $pwyKey
[[ $verbose -ge 1 ]] && echo Number of pathways to be considered: $pwyNr

pwyDBfile=$(mktemp -p $tmpdir)
echo "$pwyDB" > $pwyDBfile
Rscript $dir/prepare_batch_alignments.R $pwyDBfile $database $taxonomy $seqSrc $force_offline $update_manually $use_gene_seq $n_threads $verbose
# the final reference sequences are stored by the above R-script in "query.faa"

#----------------------#
# Calculate Alignments #
#----------------------#
touch aligner.log

if [ "$aliTool" == "blast" ]; then
    echo `blastp -version` >> aligner.log
    makeblastdb -in "$fasta" -dbtype prot -out orgdb >> aligner.log
    blastp -db orgdb -query query.faa -qcov_hsp_perc $covcutoff -num_threads $n_threads -outfmt "6 $blast_format" > alignments.tsv
fi

if [ "$aliTool" == "diamond" ]; then
    echo `diamond --version` >> aligner.log
    diamond makedb --in "$fasta" -d orgdb >> aligner.log 2>&1
    diamond blastp -d orgdb.dmnd -q query.faa \
      --threads $n_threads \
      --out alignments.tsv \
      --outfmt 6 $diamond_format \
      --query-cover $covcutoff >> aligner.log 2>&1
fi

if [ "$aliTool" == "mmseqs2" ]; then
    echo `mmseqs version` >> aligner.log
    mmseqs createdb "$fasta" targetDB >> aligner.log 2>&1
    mmseqs createdb query.faa queryDB >> aligner.log 2>&1
    mmseqs search queryDB targetDB resultDB $tmpdir --threads $n_threads -c 0.$covcutoff >> aligner.log
    mmseqs convertalis queryDB targetDB resultDB alignments.tsv \
      --format-output "$mmseqs_format" >> aligner.log 2>&1
    
    sed -Ei 's/^([^ ]+) [^\t]+/\1/' alignments.tsv
fi

# cp alignments.tsv ~/tmp/alignments.tsv # debug line

#----------------------#
# Analyse Alignments   #
#----------------------#

Rscript $dir/analyse_alignments.R $bitcutoff $identcutoff $strictCandidates $identcutoff_exception $subunit_cutoff $completenessCutoffNoHints $completenessCutoff $n_threads $vagueCutoff $verbose $pwyDBfile

#------------------------#
# Exporting result files # 
#------------------------#

# add gapseq version and sequence database status to table comments head
gapseq_version=$($dir/.././gapseq -v | head -n 1)
seqdb_version=`md5sum $dir/../dat/seq/$taxonomy/rev/sequences.tar.gz | cut -c1-7`
seqdb_date=$(stat -c %y $dir/../dat/seq/$taxonomy/rev/sequences.tar.gz | cut -c1-10)
nORFs=$(grep -c "^>" "$fasta")
nORFsMapped=$(cat nmappedORFs.tmp)
ORFcov=`echo "scale=2; $nORFsMapped*100/$nORFs" | bc`

if [ $input_mode == "nucl" ] && [ $newtranslate == "true" ]; then
    gzip -c $fasta > "$output_dir/${fastaID}.faa.gz"
    mv ${fastaID}.gff "$output_dir/${fastaID}.gff"
fi

genome_info="genome_format=${input_mode};taxonomy=${taxonomy};ORF_coverage=${ORFcov}"
[[ $taxonomy == "Bacteria" ]] && genome_info="${genome_info};gram=${gramstaining}"
[[ $input_mode == "nucl" ]] && faamd5=`md5sum $output_dir/${fastaID}.faa.gz | cut -c1-7` && genome_info="${genome_info};translation_md5=${faamd5}"
[[ $input_mode == "nucl" ]] && [[ $newtranslate == "true" ]] && genome_info="${genome_info};translation_table=${transl_table}"

sed -i "1s/^/# $gapseq_version\n/" output.tbl
sed -i "2s/^/# Sequence DB md5sum: $seqdb_version ($seqdb_date, $taxonomy)\n/" output.tbl
sed -i "3s/^/# $genome_info\n/" output.tbl
sed -i "1s/^/# $gapseq_version\n/" output_pwy.tbl
sed -i "2s/^/# Sequence DB md5sum: $seqdb_version ($seqdb_date, $taxonomy)\n/" output_pwy.tbl
sed -i "3s/^/# $genome_info\n/" output_pwy.tbl

cp aligner.log $output_dir/${fastaID}-$output_suffix-find_aligner.log
cp output.tbl $output_dir/${fastaID}-$output_suffix-Reactions.tbl
cp output_pwy.tbl $output_dir/${fastaID}-$output_suffix-Pathways.tbl




# print annotation genome coverage
#[[ $verbose -ge 1 ]] && [[ "$anno_genome_cov" = true ]] && Rscript $dir/coverage.R "$fasta" $output_dir/${fastaID}-$output_suffix-Reactions.tbl 

# cleaning
#[[ -s "$tmp_fasta" ]] && rm "$tmp_fasta"


ps -q $$ -o %cpu,%mem,args
end_time=`date +%s`
echo Running time: `expr $end_time - $start_time` s.
