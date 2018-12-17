#
# TODO: update sub2pwy to contain more linked substances
#

bitcutoff=50 # cutoff blast: min bit score
identcutoff=0   # cutoff blast: min identity
covcutoff=0 # cutoff blast: min coverage
use_alternatives=true
includeSeq=false

usage()
{
    echo "Usage"
    echo "$0 file.fasta."
    echo "  -b bit score cutoff for local alignment (default: $bitcutoff)"
    echo "  -i identity cutoff for local alignment (default: $identcutoff)"
    echo "  -c coverage cutoff for local alignment (default: $covcutoff)"
    echo "  -q Include sequences of hits in log files; default $includeSeq"
exit 1
}

OPTIND=1         # Reset in case getopts has been used previously in the shell.
while getopts "h?i:b:c:q" opt; do
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
    esac
done
shift $((OPTIND-1))
[ "$1" = "--" ] && shift
# after parsing arguments, only fasta file shoud be there
[ "$#" -ne 1 ] && { usage; }

curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
tcdb=$dir/dat/seq/tcdb.fasta
otherDB=$dir/dat/seq/transporter.fasta
subDB=$dir/dat/sub2pwy.csv
seedDB=$dir/dat/seed_transporter.tbl

# tmp working directory
fasta=$(readlink -f "$1") # save input file before changing to temporary directory
cd $(mktemp -d)

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

cat $tcdb $otherDB > all.fasta # join transporter databases
grep -e ">" all.fasta > tcdb_header
sed '1d' $subDB | awk -F ',' '{if ($8 != "NA") print $0}' > redSubDB
key=$(cat redSubDB | awk -F ',' '{if ($2 != "") print $1"|"$2; else print $1}' | paste -s -d '|') # ignore substances without linked exchange reaction
grep -wEi "$key" tcdb_header | awk '{print substr($1,2)}' > hits
#grep -wEi "arabinose" tcdb_header | awk '{print substr($1,2)}' > hits
subhits=$(grep -wEio "$key" tcdb_header | sort | uniq | paste -s -d '|')
allDBsubs=$(cat redSubDB | grep -wEi "$subhits" | awk -F ',' '{if ($2 != "") print $2; else print $1}' | sort | paste -s -d ';') # list of substances that are covered by DB
#apt install exonerate
fastaindex all.fasta tcdb.idx 
fastafetch -f all.fasta -i tcdb.idx -Fq <(sort -u hits ) > tcdbsmall.fasta

#
makeblastdb -in $fasta -dbtype nucl -out orgdb >/dev/null

#tblastn -db orgdb -query tcdbsmall.fasta -outfmt "6 $blast_format" > out 
cat tcdbsmall.fasta | parallel --will-cite --block 500k --recstart '>' --pipe tblastn -db orgdb -qcov_hsp_perc $covcutoff -outfmt \'"6 $blast_format"\' -query - > out

IDtcdb=$(cat out | awk -v bitcutoff=$bitcutoff -v identcutoff=$identcutoff -v covcutoff=$covcutoff '{if ($2>=identcutoff && $5>=covcutoff && $4>=bitcutoff) print $1}' | cut -d "|" -f 3 | sort | uniq) 

TC[1]="1.Channels and pores"
TC[2]="2.Electrochemical potential-driven transporters"
TC[3]="3.Primary active transporters"
TC[4]="4.Group translocators"


for id in $IDtcdb
do
    descr=$(grep $id tcdb_header | grep -P "(?<= ).*" -o) # extract description
    tc=$(grep $id tcdb_header | grep -Pw "([1-4]\\.[A-Z]+\\.[0-9]+\\.[0-9]+)" -o) # ATTENTION: only TC numbers starting with 1,2,3,4 are selected (others are electron carrier and accessoirs)
    i=$(echo $tc | head -c1) # get first character of TC number => transporter type
    type=${TC[$i]}
    [ -z "$tc" ] && continue
    subl=$(echo "$descr" | grep -wEio "$key" | grep -if - redSubDB | awk -F ',' '{ if ($2 != "") print $2; else print $1}') # could be more then one hit (e.g. transporter with two substances)
    for subst in $subl
    do
        #echo -e $subst"\t"$type"\t"$tc"\t"$descr >> newTransporter.tbl
        exmetall=$(cat $subDB | grep -w "$subst" | awk -F ',' '{if ($8 != "NA") print $8}')
        for exmet in $exmetall
        do
            exid=$(cat $subDB | grep -w "$subst" | awk -F ',' '{if ($7 != "NA") print $7}')
            sublist="$sublist;$subst"
            echo $id,$tc,$subst,$exmet,$exid,$type,$descr >> `echo "$subst" | tr '[:upper:]' '[:lower:]'`
            cat $seedDB | awk -F '\t' -v type="$type" -v exmet="$exmet" '{if($3==type && $4==exmet) print $0}' > tr_cand
            if [ -s tr_cand ]; then
                rea=$(cat tr_cand | awk -F '\t' '{print $1}')
                #echo -e "\t"Found: $rea
                cand="$cand $rea $exid"
                sublist2="$sublist2;$subst"
                rea_export=$(echo "$rea" | tr '\n' ',' | sed 's/,$//g')
                cat out | grep -w $id | sort -rgk 4,4 | head -1 | awk -v id="$id" -v tc="$tc" -v subst="$subst" -v exid="$exid" -v rea=$rea_export '{print id"\t"tc"\t"subst"\t"exid"\t"rea"\t"$0}' >> transporter.tbl 
            fi
        done
    done
done

echo -e "\nFound transporter for:"
echo ${sublist:1} | tr '[:upper:]' '[:lower:]' | tr ';' '\n' | sort | uniq -c
echo ${sublist:1} | tr '[:upper:]' '[:lower:]' | tr ';' '\n' | sort | uniq  > hit1

echo -e "\nFound transporter and import reactions for:"
echo ${sublist2:1} | tr '[:upper:]' '[:lower:]' | tr ';' '\n' | sort | uniq > hit2
cat hit2

echo -e "\nNo transport reactions found in database for:"
for sub in $(comm -23 hit1 hit2)
do
    cat "$sub" | awk -F ',' '{print $3, $6}' | sort | uniq
    for exmet in $(cat "$sub" | cut -d ',' -f 4 | sort | uniq)
    do
        alter=$(cat $seedDB | awk -F '\t' -v exmet="$exmet" '{if($4==exmet) print $1}' | sort | uniq)
        if [ -n "$alter" ]; then
            echo -e "\tAlternatives:" $alter
            if [ "$use_alternatives" = true ] ; then
                cand="$cand $alter $exid"
            fi
        fi
    done
done

echo -e "\nNo transporter found for compounds (existing transporter/exchanges should be removed?):"
echo $allDBsubs | tr '[:upper:]' '[:lower:]' | tr ';' '\n' | sort > hit3
cat hit1 | tr '\n' '|' | rev | cut -c 2- | rev | grep -f - -iE $subDB | cut -d ',' -f1,2 | tr '[:upper:]' '[:lower:]' | tr ',' '\n' | sort > hit4 # expand hits with substance alternative names to avoid false positive reporting
comm -13 hit4 hit3

cand="$(echo $cand | tr ' ' '\n' | sort | uniq | tr '\n' ' ')"
#echo -e "\nReactions to be added:"
#echo $cand
echo $cand > newTransporter.lst
cp newTransporter.lst $curdir/${fastaid}-Transporter.lst
#cp newTransporter.tbl $curdir/${fastaid}-Transporter.tbl
cp transporter.tbl $curdir/${fastaid}-Transporter.tbl
