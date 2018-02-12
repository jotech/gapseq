#
# TODO: update sub2pwy to contain more linked substances
#

curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
#fasta=/home/jo/uni/gapseq/dat/myb71.fna
fasta=$(readlink -f $1)
tmpvar=$(basename $fasta)
fastaid=${tmpvar%.*}
tcdb=/home/jo/uni/gapseq/dat/tcdb.fasta
subDB=/home/jo/uni/gapseq/dat/sub2pwy.csv
seedDB=/home/jo/uni/gapseq/dat/reftrdb_seed.csv


# tmp working directory
cd $(mktemp -d)


grep -e ">" $tcdb > tcdb_header
sed '1d' $subDB | awk -F ',' '{if ($8 != "NA") print $0}' > redSubDB
key=$(cat redSubDB | awk -F ',' '{if ($2 != "") print $1"|"$2; else print $1}' | paste -s -d '|') # ignore substances without linked exchange reaction
allsubs=$(cat redSubDB | awk -F ',' '{if ($2 != "") print $2; else print $1}' | paste -s -d ';') 
grep -wEi "$key" tcdb_header | awk '{print substr($1,2)}' > hits
subhits=$(grep -wEio "$key" tcdb_header | sort | uniq | paste -s -d '|')
allDBsubs=$(cat redSubDB | grep -wEi "$subhits" | awk -F ',' '{if ($2 != "") print $2; else print $1}' | sort | paste -s -d ';') # list of substances that are covered by DB
#apt install exonerate
fastaindex $tcdb tcdb.idx 
fastafetch -f $tcdb -i tcdb.idx -Fq <(sort -u hits ) > tcdbsmall.fasta

#
makeblastdb -in $fasta -dbtype nucl -out orgdb >/dev/null

tblastn -db orgdb -query tcdbsmall.fasta -outfmt '6 qseqid sseqid pident evalue bitscore stitle' > out 

IDtcdb=$(cat out | awk '{if ($5>=50 && $3>=50) print $1}' | awk -F "|" '{print $3}') 

TC[1]="1.Channels and pores"
TC[2]="2.Electrochemical potential-driven transporters"
TC[3]="3.Primary active transporters"
TC[4]="4.Group translocators"


for id in $IDtcdb
do
    descr=$(grep $id tcdb_header | grep -P "(?<= ).*(?=OS)" -o) # extract description
    tc=$(grep $id tcdb_header | grep -Pw "([0-9]+.[A-Z]+.[0-9]+.[0-9]+)" -o)
    [ -z "$descr" ] && descr=$(grep $id tcdb_header | grep -P "(?<= ).*(?= - )" -o)
    subl=$(echo "$descr" | grep -wEio "$key" | grep -if - redSubDB | awk -F ',' '{ if ($2 != "") print $2; else print $1}') # could be more then one hit (e.g. transporter with two substances)
    for sub in $subl
    do
        exmet=$(cat $subDB | grep -w "$sub" | awk -F ',' '{if ($8 != "NA") print $8}')
        exid=$(cat $subDB | grep -w "$sub" | awk -F ',' '{if ($7 != "NA") print $7}')
        sublist="$sublist;$sub"
        echo $id $tc $sub $exmet $exid $descr
        i=$(echo $tc | head -c1) # get first character of TC number => transporter type
        type=${TC[$i]}
        cat $seedDB | awk -F '\t' -v type="$type" -v exmet="$exmet" '{if($8==type && $3==exmet) print $0}' > tr_cand
        if [ -s tr_cand ]; then
            rea=$(cat tr_cand | awk -F '\t' '{print $1}')
            echo -e "\t"Found: $rea
            cand="$cand $rea $exid"
            sublist2="$sublist2;$sub"
        fi
    done
done

echo -e "\nFound transporter for:"
echo ${sublist:1} | tr '[:upper:]' '[:lower:]' | tr ';' '\n' | sort | uniq -c
echo ${sublist:1} | tr '[:upper:]' '[:lower:]' | tr ';' '\n' | sort | uniq  > hit1

echo -e "\nFound transporter and import reactions for:"
echo ${sublist2:1} | tr '[:upper:]' '[:lower:]' | tr ';' '\n' | sort | uniq > hit2
cat hit2

echo -e "\nDo not found transport reactions in database for:"
comm -23 hit1 hit2

echo -e "\n No transporter found for compounds (existing transporter/exchanges should be removed?):"
echo $allDBsubs | tr ';' '\n' | sort > hit3
comm -13 hit1 hit3

cand="$(echo $cand | tr ' ' '\n' | sort | uniq | tr '\n' ' ')"
echo -e "\nReactions to be added:"
echo $cand
echo $cand > newTransporter.lst
cp newTransporter.lst $curdir/${fastaid}Transporter.lst