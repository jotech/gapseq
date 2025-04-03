#!/bin/bash

start_time=`date +%s`
bitcutoff=50
identcutoff=0   # cutoff blast: min identity
covcutoff=75 # cutoff blast: min coverage
use_alternatives=true
includeSeq=false
use_parallel=true
only_met=""
verbose=1
user_temp=false
input_mode="auto"
output_dir=.
OS=$(uname -s)
if [ "$OS" = "Darwin" -o "$OS" = "FreeBSD" ]; then
	n_threads=$(sysctl hw.ncpu|cut -f2 -d' ')
else
	n_threads=`grep -c ^processor /proc/cpuinfo`
fi

usage()
{
    echo "Usage"
    echo "$0 file.fasta."
    echo "  -b bit score cutoff for local alignment (default: $bitcutoff)"
    echo "  -i identity cutoff for local alignment (default: $identcutoff)"
    echo "  -c coverage cutoff for local alignment (default: $covcutoff)"
    echo "  -q Include sequences of hits in log files; default $includeSeq"
    echo "  -k do not use parallel"
    echo "  -m only check for this keyword/metabolite (default: all)"
    echo "  -f Path to directory, where output files will be saved (default: current directory)"
    echo "  -v Verbose level, 0 for nothing, 1 for full (default $verbose)"
    echo "  -T Set user-defined temporary folder (default: $user_temp)"
    echo "  -M Input genome mode. Either 'nucl' or 'prot' (default '$input_mode')"
    echo "  -K Number of threads for sequence alignments. If option is not provided, number of available CPUs will be automatically determined."
    echo ""
exit 1
}

OPTIND=1         # Reset in case getopts has been used previously in the shell.
while getopts "h?i:b:c:qkm:f::v:T:M:K:" opt; do
    case "$opt" in
    h|\?)
        usage
        exit 0
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
    q)
        includeSeq=true
        ;;
    k)
        use_parallel=false
        ;;
    m)
        only_met=$OPTARG
        ;;
    f)
        output_dir=$OPTARG
        ;;
    v)  
        verbose=$OPTARG
        ;;
    T)
        user_temp=true
        user_temp_folder=$OPTARG
        ;;
    M)
        input_mode=$OPTARG
        ;;
    K)
        n_threads=$OPTARG
        ;;
    esac
done
shift $((OPTIND-1))
[ "$1" = "--" ] && shift
# after parsing arguments, only fasta file shoud be there
[ "$#" -ne 1 ] && { usage; }

curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
tcdb=$dir/../dat/seq/tcdb.fasta
otherDB=$dir/../dat/seq/transporter.fasta
subDB=$dir/../dat/subex.tbl
seedDB=$dir/../dat/seed_transporter.tbl
customDB=$dir/../dat/seed_transporter_custom.tbl
tcdb_sub=$dir/../dat/tcdb_substrates.tbl
tcdb_custom=$dir/../dat/tcdb_custom.tbl

# set output directory
case $output_dir in
    /*)
        # echo "absolute path"
        output_dir=$output_dir
        ;;
    ~*)
        # echo "relative to $HOME directory"
        output_dir="${output_dir/#\~/$HOME}"
        ;;
    *)
        # echo "relative path to current directory"
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
echo $tmpdir
cd $tmpdir

cat $seedDB $customDB > allDB

# Get fasta file
if [[ $fasta == *.gz ]]; then # in case fasta is in an archive
    tmp_fasta=$(basename "${fasta}" .gz)
    gunzip -c $fasta > $tmp_fasta
    fasta=$tmp_fasta
fi
[[ ! -s $fasta ]] && { echo Invalid file: $1; exit 0; }

# Determine if fasta is nucl or prot
if [ $input_mode == "auto" ]; then
    input_mode=`$dir/./nuclprot.sh $fasta`
    
    if [ $input_mode == "prot" ]; then
        echo "Protein fasta detected."
    else
        echo "Nucleotide fasta detected."
    fi
    
fi

tmpvar=$(basename $fasta)
fastaid=${tmpvar%.*}

# blast format
if [ "$includeSeq" = true ]; then
    blast_format="qseqid pident evalue bitscore qcovs stitle sstart send sseq"
else
    blast_format="qseqid pident evalue bitscore qcovs stitle sstart send"
fi

cat $tcdb $otherDB > all.fasta # join transporter databases
grep -e ">" all.fasta | sort > fasta_header
sed '1d' $subDB | awk -F '\t' '{if ($4 != "") print $0}' > redSubDB
cat $tcdb_sub $tcdb_custom > tcdb_all

# custom fasta_header corrections to avoid mismatches
sed -i 's/amino butyrate/aminobutyrate/g' fasta_header
sed -i 's/GABA butyrate/GABA/g' fasta_header

if [[ -n "$only_met" ]]; then
    cat redSubDB | grep -wi "$only_met" | awk -F '\t' '{n=split($1,a,";"); for(i=0;++i<=n;){print a[i]"\t"$4}}' > SUBid
    echo Search keys: `cat SUBid | cut -f1 | paste -s -d ';'`
else
    cat redSubDB | awk -F '\t' '{n=split($1,a,";"); for(i=0;++i<=n;){print a[i]"\t"$4}}' > SUBid
fi
cat SUBid | cut -f1 | sort | uniq > SUBkey
grep -Fiwf SUBkey tcdb_all | cut -f1 > TCkey
grep -Fivf TCkey fasta_header > fasta_header.noTCkey
grep -Fivf SUBkey fasta_header.noTCkey > fasta_header.noKey
comm -23 fasta_header fasta_header.noKey > fasta_header.small
awk 'BEGIN{while((getline<"fasta_header.small")>0)l[$1]=1}/^>/{f=l[$1]}f' all.fasta > small.fasta 


makeblastdb -in $fasta -dbtype $input_mode -out orgdb >/dev/null
if [ "$input_mode" == "nucl" ]; then
    tblastn -db orgdb -qcov_hsp_perc $covcutoff -outfmt "6 $blast_format" -query small.fasta > out
fi
if [ "$input_mode" == "prot" ]; then
    blastp -db orgdb -qcov_hsp_perc $covcutoff -num_threads $n_threads -outfmt "6 $blast_format" -query small.fasta > out
fi



cat out | awk -v identcutoff=$identcutoff -v covcutoff=$covcutoff '{if ($2>=identcutoff && $5>=covcutoff) print $0}'> blasthits
TC_blasthits=$(grep -Pwo "([1-4]\\.[A-z]\\.[0-9]+\\.[0-9]+\\.[0-9]+)" blasthits | sort | uniq)
TC_types[1]="1.Channels and pores"
TC_types[2]="2.Electrochemical potential-driven transporters"
TC_types[3]="3.Primary active transporters"
TC_types[4]="4.Group translocators"

touch DBfound DBfound_bitscore DBmissing DBmissing_bitscore transporter.tbl transporter_bitscore.tbl
IFS=$'\n' # only splits by newline in for-loop
for tc in $TC_blasthits
do
    i=$(echo $tc | head -c1) # get first character of TC number => transporter type
    type=${TC_types[$i]}
    substrates=$(grep -F $tc tcdb_all | cut -f2)
    sub_hit=$(echo "$substrates" | grep -woiFf SUBkey | sort | uniq)
    if [[ -z "$sub_hit" ]]; then
        sub_hit=$(grep -F $tc fasta_header.small | grep -woiFf SUBkey | sort | uniq)
    fi
    if [[ -z "$sub_hit" ]]; then
        echo -e "$substrates\t$tc\t$type" >> SUBmissing
        continue
    fi
    best_hit=$(grep -F $tc blasthits | sort -rgk 4,4 | head -1)
    hit_good=$(echo "`echo "$best_hit" | cut -f4` >= $bitcutoff" | bc)

    for sub in $sub_hit
    do 
        exid_all=$(awk -F '\t' -v subl="$sub" 'tolower($1)==tolower(subl) {print $2}' SUBid | sort | uniq) # get id of substance to avoid mismatches
        for exid in $exid_all
        do 
            exmet=$(echo "$exid" | sed 's/EX_//g' | sed 's/_e0/\[e0\]/g')
            candidates=$(cat allDB | awk -F '\t' -v type="$type" -v exmet="$exmet" '{if($3==type && $4==exmet) print $0}')
            if [[ -n "$candidates" ]]; then
                rea=$(echo "$candidates" | cut -f1)
                rea_export=$(echo "$rea" | tr '\n' ',' | sed 's/,$//g')
                exid_export=$(echo "$exid" | tr '\n' ',' | sed 's/,$//g')
                echo -e "`echo $best_hit | cut -f1`\t$tc\t$sub\t$exid_export\t$rea_export\t$best_hit\ttransporter" >> transporter.tbl
                echo -e "$sub\t$exmet\t$tc\t$type" >> DBfound
                if [[ $hit_good -eq 1 ]]; then
                    [[ verbose -ge 1 ]] && echo $sub_hit $tc $type 
                    [[ -n "$only_met" ]] && { echo -e "\t"Found reaction: $rea; } 
                    echo -e "`echo $best_hit | cut -f1`\t$tc\t$sub\t$exid_export\t$rea_export\t$hit_good\ttransporter" >> transporter_bitscore.tbl
                    echo -e "$sub\t$exmet\t$tc\t$type" >> DBfound_bitscore
                fi
            else
                echo -e "$sub\t$exmet\t$tc\t$type" >> DBmissing
                [[ $hit_good -eq 1 ]] && echo -e "$sub\t$exmet\t$tc\t$type" >> DBmissing_bitscore
            fi
        done
    done
done

cat DBfound | cut -f1 | sort | uniq -i > IDfound
cat DBfound_bitscore | cut -f1 | sort | uniq -i > IDfound_bitscore
[[ verbose -ge 1 ]] && { echo -e "\nFound transporter for substances:"; cat IDfound_bitscore; }
cat DBmissing | cut -f1 | sort | uniq -i > IDmissing
cat DBmissing_bitscore | cut -f1 | sort | uniq -i > IDmissing_bitscore
[[ verbose -ge 1 ]] && { echo -e "\nFound transporter without reaction in DB:"; comm -23 IDmissing_bitscore IDfound_bitscore; }


if [ "$use_alternatives" = true ] ; then
    [[ verbose -ge 1 && -s IDmissing ]] && echo -e "\nAdd alternative transport reactions:"
    for missing in `comm -23 IDmissing IDfound`
    do
        [[ verbose -ge 1 ]] && echo $missing
        exid_all=$(awk -F '\t' -v subl="$missing" 'tolower($1)==tolower(subl) {print $2}' SUBid | sort | uniq)
        for exid in $exid_all
        do
            exmet=$(echo "$exid" | sed 's/EX_//g' | sed 's/_e0/\[e0\]/g')
            candidates=$(cat allDB | awk -F '\t' -v exmet="$exmet" '{if($4==exmet) print $0}')
            rea=$(echo "$candidates" | cut -f1)
            rea_export=$(echo "$rea" | tr '\n' ',' | sed 's/,$//g')
            exid_export=$(echo "$exid" | tr '\n' ',' | sed 's/,$//g')
            tc=$(grep $missing DBmissing | cut -f3 | paste -s -d '|')
            best_alter=$(grep -E $tc blasthits | sort -rgk 4,4 | head -1)
            echo -e "`echo $best_alter | cut -f1`\t`echo $best_alter|grep -Pwo "([1-4]\\.[A-z]\\.[0-9]+\\.[0-9]+\\.[0-9]+)"`\t$missing\t$exid_export\t$rea_export\t$best_alter\talt-transporter" >> transporter.tbl
        done
    done
fi
unset IFS
echo Total number of found substance transporter for $fastaid: `cat transporter_bitscore.tbl | cut -f4 | sort | uniq | wc -l`

cp transporter.tbl $output_dir/${fastaid}-Transporter.tbl
[[ -s transporter.tbl ]] && echo "id tc sub exid rea $blast_format comment" | tr ' ' '\t' | cat - transporter.tbl | awk '!a[$0]++' > $output_dir/${fastaid}-Transporter.tbl # add header and remove duplicates

# add gapseq vesion and sequence database status to table comments head
gapseq_version=$($dir/.././gapseq -v | head -n 1)
seqdb_version=$(md5sum $dir/../dat/seq/transporter.fasta | cut -c1-7)
seqdb_date=$(stat -c %y $dir/../dat/seq/transporter.fasta | cut -c1-10)
sed -i "1s/^/# $gapseq_version\n/" $output_dir/${fastaid}-Transporter.tbl
sed -i "2s/^/# Transporter sequence DB md5sum: $seqdb_version ($seqdb_date)\n/" $output_dir/${fastaid}-Transporter.tbl

# finishing
end_time=`date +%s`
echo Running time: `expr $end_time - $start_time` s.
