#!/bin/bash

#[ $# -ne 2 ] && { echo "Usage: $0 file.fasta model.sbml"; exit 1; }

# paths and variables
#fasta=$2
fasta=~/uni/dat/seq/Celegans_Microbiome/genomes/all/MYb71.fasta
csource="Galactose"
path=$(readlink -f "$0")
dir=$(dirname "$path")
pwydb=$dir/dat/pwydf.csv
seqpath=$dir/dat/seq/
bitcutoff=50 # cutoff blast: min bit score
readb=$dir/dat/vmh_reactions.csv


# tmp working directory
#cd $(mktemp -d)
cd /tmp/tmp.OKh8bqa0QC

# create blast database
makeblastdb -in $fasta -dbtype nucl -out orgdb >/dev/null

hit=$(grep $csource $pwydb)

echo Check for carbon source: $csource
for i in `seq 1 $(echo $hit | wc -l)`
do
    count=0
    countex=0
    countdb=0
    pwy=$(echo $hit | awk -v i=$i 'NR==i {print $3}')
    ecs=$(echo $hit | awk -v i=$i 'NR==i {print $4}')
    reaids=$(echo $hit | awk -v i=$i 'NR==i {print $5}')
    echo Checking for pathway $pwy
    for ec in $(echo $ecs | tr "," "\n")
    do 
        ((count++))
        query=$seqpath$ec.fasta
        out=$pwy-$ec.blast
        rea=$(echo $reaids | tr "," "\n" | awk -v j=$count 'NR==j {print $1}')
        tblastn -db orgdb -query $query -outfmt '6 qseqid sseqid pident evalue bitscore stitle' >$out 
        if [ -s $out ]; then
            hits=$(cat $out | awk -v bitcutoff=$bitcutoff '$5>=bitcutoff {print $1}' | wc -l)
            if [ $hits -gt 0 ]; then
                echo Found reaction: $rea $ec
                dbhit=$(grep -w $ec $readb | awk -F ',' '{print $1}')
                if [ -n "$dbhit" ]; then
                    echo Candidate reaction found for import: $dbhit
                    ((countdb++))
                else
                    echo No candidate reaction found for import: $rea
                fi
                ((countex++))
            else
                echo No significant blast hits found for reaction: $rea 
            fi
        else
            echo No blast hits found for reaction: $rea 
        fi
    done
    echo Pathway completness: $countex/$count
    echo Candidate reactions found in db: $countdb
done

