#
# TODO: update sub2pwy to contain more linked substances
#

start_time=`date +%s`

bitcutoff=50 # cutoff blast: min bit score
identcutoff=0   # cutoff blast: min identity
covcutoff=75 # cutoff blast: min coverage
use_alternatives=true
includeSeq=false
use_parallel=true
only_met=""
verbose=1

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
    echo "  -v Verbose level, 0 for nothing, 1 for full (default $verbose)"
exit 1
}

OPTIND=1         # Reset in case getopts has been used previously in the shell.
while getopts "h?i:b:c:qkm:v:" opt; do
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
    v)  
        verbose=$OPTARG
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
subDB=$dir/../dat/sub2pwy.csv
seedDB=$dir/../dat/seed_transporter.tbl
customDB=$dir/../dat/seed_transporter_custom.tbl

# tmp working directory
fasta=$(readlink -f "$1") # save input file before changing to temporary directory
cd $(mktemp -d)

cat $seedDB $customDB > allDB

# Get fasta file
if [[ $fasta == *.gz ]]; then # in case fasta is in a archive
    tmp_fasta=$(basename "${fasta}" .gz)
    gunzip -c $fasta > $tmp_fasta
    fasta=$tmp_fasta
fi
[[ ! -s $fasta ]] && { echo Invalid file: $1; exit 0; }

tmpvar=$(basename $fasta)
fastaid=${tmpvar%.*}

# blast format
if [ "$includeSeq" = true ]; then
    blast_format="qseqid pident evalue bitscore qcovs stitle sstart send sseq"
else
    blast_format="qseqid pident evalue bitscore qcovs stitle sstart send"
fi

cat $tcdb $otherDB > all.fasta.tmp # join transporter databases
perl -ne 'if (/>(.*?)\s+(.*)/){push(@{$hash{$1}},$2) ;}}{open(I, "<","all.fasta.tmp");while(<I>){if(/>(.*?)\s+/){ $t = 0; next if $h{$1}; $h{$1} = 1 if $hash{$1}; $t = 1; chomp; print $_ . " @{$hash{$1}}\n"}elsif($t==1){print $_} } close I;' all.fasta.tmp > all.fasta # remove duplicate entries by ID
sed -i "s/\(>.*\)/\L\1/" all.fasta # header to lower case
grep -e ">" all.fasta > tcdb_header
sed '1d' $subDB | awk -F '\t' '{if ($8 != "NA") print $0}' > redSubDB
[[ -n "$only_met" ]] && { cat redSubDB | grep -i $only_met > redSubDB.tmp; mv redSubDB.tmp redSubDB; }
[[ ! -s redSubDB ]] && { echo keyword/metabolite not found; exit 1; }
key=$(cat redSubDB | awk -F '\t' '{if ($2 != "") print $1"\n"$2; else print $1}' | sort | uniq | paste -s -d '|') # ignore substances without linked exchange reaction
[[ -n "$only_met" ]] && { echo Search keys: "$key"; } 
grep -wEi "$key" tcdb_header | awk '{print substr($1,2)}' > hits
#grep -wEio "$key" tcdb_header | sort | uniq -c # status print
[[ -n "$only_met" ]] && { echo -e "\n"Found transporter sequences:; cat hits; } 
subhits=$(grep -wEio "$key" tcdb_header | sort | uniq | paste -s -d '|')
allDBsubs=$(cat redSubDB | grep -wEi "$subhits" | awk -F '\t' '{if ($2 != "") print $2; else print $1}' | sort | uniq | paste -s -d ';') # list of substances that are covered by DB
[[ -n "$only_met" ]] && { echo -e "\n"Covered metabolites: $allDBsubs; } 
touch transporter.tbl
#apt install exonerate
fastaindex all.fasta tcdb.idx 
fastafetch -f all.fasta -i tcdb.idx -Fq <(sort -u hits ) > tcdbsmall.fasta
#
makeblastdb -in $fasta -dbtype nucl -out orgdb >/dev/null

#tblastn -db orgdb -query tcdbsmall.fasta -outfmt "6 $blast_format" > out 
if ! [ -x "$(command -v parallel)" ] || [ "$use_parallel" = false ]; then # try to use parallelized version
    tblastn -db orgdb -qcov_hsp_perc $covcutoff -outfmt "6 $blast_format" -query tcdbsmall.fasta > out
else
    cat tcdbsmall.fasta | parallel --gnu --will-cite --block 500k --recstart '>' --pipe tblastn -db orgdb -qcov_hsp_perc $covcutoff -outfmt \'"6 $blast_format"\' -query - > out
fi

IDtcdb=$(cat out | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff -v covcutoff=$covcutoff '{if ($2>=identcutoff && $5>=covcutoff) print $1}' | cut -d "|" -f 3 | sort | uniq) 

TC[1]="1.Channels and pores"
TC[2]="2.Electrochemical potential-driven transporters"
TC[3]="3.Primary active transporters"
TC[4]="4.Group translocators"

[[ -n "$only_met" ]] && { echo -e "\n"Blast hits:; } 
for id in $IDtcdb
do
    descr=$(grep $id tcdb_header | grep -P "(?<= ).*" -o) # extract description
    [[ -n "$only_met" ]] && { echo -e "\t\n"$descr; } 
    tc=$(grep $id tcdb_header | grep -Pw "([1-4]\\.[A-z]\\.[0-9]+)" -o | uniq | tr '\n' ',' | sed 's/,$//g') # ATTENTION: only TC numbers starting with 1,2,3,4 are selected (others are electron carrier and accessoirs)
    i=$(echo $tc | head -c1) # get first character of TC number => transporter type
    type=${TC[$i]}
    [ -z "$tc" ] && continue
    subl=$(echo "$descr" | grep -wEio "$key" | tr ' ' '_' | sort | uniq)
    [[ -n "$only_met" ]] && { echo -e "\t"$id $tc $subl; } 
    #echo $id $descr $tc $i $type $subl
    for subst in $subl
    do
        subst=$(echo "$subst" | tr '_' ' ') # get space character back (was changed for look)
        exmetall=$(cat $subDB | awk -F '\t' -v subst="$subst" '{if ( (tolower($1) == tolower(subst) || tolower($2) == tolower(subst)) && $8 != "NA" ) print $8}')
        for exmet in $exmetall
        do
            exid=$(cat $subDB | grep -w "$subst" | awk -F '\t' '{if ($7 != "NA") print $7}')
            sublist="$sublist;$subst"
            echo $id,$tc,$subst,$exmet,$exid,$type,$descr >> "`echo "$subst" | tr '[:upper:]' '[:lower:]'`"
            cat allDB | awk -F '\t' -v type="$type" -v exmet="$exmet" '{if($3==type && $4==exmet) print $0}' > tr_cand
            if [ -s tr_cand ]; then
                rea=$(cat tr_cand | awk -F '\t' '{print $1}')
                [[ -n "$only_met" ]] && { echo -e "\t"Found reaction: $rea; } 
                #echo -e "\t"Found: $rea
                cand="$cand $rea $exid"
                sublist2="$sublist2;$subst"
                rea_export=$(echo "$rea" | tr '\n' ',' | sed 's/,$//g')
                exid_export=$(echo "$exid" | tr '\n' ',' | sed 's/,$//g')
                cat out | grep -w $id | sort -rgk 4,4 | head -1 | awk -v id="$id" -v tc="$tc" -v subst="$subst" -v exid="$exid_export" -v rea=$rea_export '{print id"\t"tc"\t"subst"\t"exid"\t"rea"\t"$0}' >> transporter.tbl 
            fi
        done
    done
done

[[ verbose -ge 1 ]] && echo -e "\nFound transporter for:"
[[ verbose -ge 1 ]] && echo ${sublist:1} | tr '[:upper:]' '[:lower:]' | tr ';' '\n' | sort | uniq -c
echo ${sublist:1} | tr '[:upper:]' '[:lower:]' | tr ';' '\n' | sort | uniq  > hit1

[[ verbose -ge 1 ]] && echo -e "\nFound transporter and import reactions for:"
echo ${sublist2:1} | tr '[:upper:]' '[:lower:]' | tr ';' '\n' | sort | uniq > hit2
[[ verbose -ge 1 ]] && cat hit2

[[ verbose -ge 1 ]] && echo -e "\nNo transport reactions found in database for:"
for sub in $(comm -23 hit1 hit2 | tr ' ' '_')
do
    sub=$(echo "$sub" | tr '_' ' ')
    [[ verbose -ge 1 ]] && { cat "$sub" | awk -F ',' '{print $3, $6}' | sort | uniq; }
    id=$(cat "$sub" | awk -F ',' '{print $1}' | sort | uniq | tr '\n' '|' | rev | cut -c2- | rev)
    tc=$(cat "$sub" | awk -F ',' '{print $6}' | sort | uniq)
    for exmet in $(cat "$sub" | cut -d ',' -f 4 | sort | uniq)
    do
        alter=$(cat allDB | awk -F '\t' -v exmet="$exmet" '{if($4==exmet) print $1}' | sort | uniq)
        if [ -n "$alter" ]; then
            [[ verbose -ge 1 ]] && echo -e "\tAlternatives:" $alter
            if [ "$use_alternatives" = true ] ; then
                cand="$cand $alter $exid"
                alter_str=$(echo "$alter" | tr '\n' ',' | rev | cut -c2- | rev) 
                cat out | grep -wE $id | sort -rgk 4,4 | head -1 | awk -v id="$id" -v tc="alternative" -v subst="$sub" -v exid="$exmet" -v rea="$alter_str" '{print id"\t"tc"\t"subst"\t"exid"\t"rea"\t"$0}' >> transporter.tbl 
            fi
        fi
    done
done

[[ verbose -ge 1 ]] && echo -e "\nNo transporter found for compounds (existing transporter/exchanges should be removed?):"
echo $allDBsubs | tr '[:upper:]' '[:lower:]' | tr ';' '\n' | sort > hit3
cat hit1 | tr '\n' '|' | rev | cut -c 2- | rev | grep -f - -iE $subDB | cut -d '	' -f1,2 | tr '[:upper:]' '[:lower:]' | tr '\t' '\n' | tr ',' '\n' | sort | uniq > hit4 # expand hits with substance alternative names to avoid false positive reporting
[[ verbose -ge 1 ]] && comm -13 hit4 hit3

cand="$(echo $cand | tr ' ' '\n' | sort | uniq | tr '\n' ' ')"
#echo -e "\nReactions to be added:"
#echo $cand
echo $cand > newTransporter.lst
#cp newTransporter.lst $curdir/${fastaid}-Transporter.lst
cp transporter.tbl $curdir/${fastaid}-Transporter.tbl
[[ -s transporter.tbl ]] && echo "id tc sub exid rea $blast_format" | tr ' ' '\t' | cat - transporter.tbl | awk '!a[$0]++' > $curdir/${fastaid}-Transporter.tbl # add header and remove duplicates

end_time=`date +%s`
echo Running time: `expr $end_time - $start_time` s.
