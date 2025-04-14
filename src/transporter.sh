#!/bin/bash

start_time=`date +%s`
bitcutoff=50
identcutoff=0   # cutoff blast: min identity
covcutoff=75 # cutoff blast: min coverage
nouse_alternatives=false
includeSeq=false
only_met=""
verbose=1
user_temp=false
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
    echo "gapseq - transporter prediction"
    echo ""
    echo "Usage"
    echo "$0 file.fasta"
    echo "  -b bit score cutoff for local alignment (default: $bitcutoff)"
    echo "  -i identity cutoff for local alignment (default: $identcutoff)"
    echo "  -c coverage cutoff for local alignment (default: $covcutoff)"
    echo "  -q Include sequences of hits in results table"
    echo "  -a Do not use alternative transport reaction types if no database reaction exists for the specific TC type of the reference sequence"
    echo "  -k Do not use parallel (Deprecated: use '-K 1' instead to disable multi-threading.)"
    echo "  -m only check for this keyword/metabolite (default: all)"
    echo "  -f Path to directory, where output files will be saved (default: current directory)"
    echo "  -v Verbose level, 0 for nothing, 1 for full (default $verbose)"
    echo "  -T Set user-defined temporary folder (default: $user_temp)"
    echo "  -M Input genome mode. Either 'nucl' or 'prot' (default '$input_mode')"
    echo "  -K Number of threads for sequence alignments. If option is not provided, number of available CPUs will be automatically determined."
    echo "  -A Tool to be used for sequence alignments (blast, mmseqs2, diamond; default: $aliTool)"
    echo ""
exit 1
}

OPTIND=1         # Reset in case getopts has been used previously in the shell.
while getopts "h?i:b:c:qakm:f::v:T:M:K:A:" opt; do
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
        n_threads=1
        echo "DEPRECATION NOTICE: Option '-k' is deprecated. To disable multi-threading use '-K 1' instead."
        ;;
    a)
        nouse_alternatives=true
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
[[ $verbose -ge 1 ]] && echo $tmpdir
cd $tmpdir

cat $seedDB $customDB > allDB

# Get fasta file
if [[ $fasta == *.gz ]]; then # in case fasta is in an archive
    tmp_fasta=$(basename "${fasta}" .gz)
    gunzip -c $fasta > $tmp_fasta
    fasta=$tmp_fasta
fi
[[ ! -s $fasta ]] && { echo Invalid file: $1; exit 0; }
tmpvar=$(basename "$fasta")
fastaID="${tmpvar%.*}"

# Determine if fasta is nucl or prot
if [ $input_mode == "auto" ]; then
    input_mode=`$dir/./nuclprot.sh $fasta`
    
    if [ $input_mode == "prot" ]; then
        echo "Protein fasta detected."
    else
        echo "Nucleotide fasta detected."
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

# alignment statistics format
if [ "$includeSeq" = true ]; then
    blast_format="qseqid pident evalue bitscore qcovs stitle sstart send sseq"
    mmseqs_format="qheader,pident,evalue,bits,qcov,theader,tstart,tend,qseq"
    diamond_format="qseqid pident evalue bitscore qcovhsp stitle sstart send sseq"
else
    blast_format="qseqid pident evalue bitscore qcovs stitle sstart send"
    mmseqs_format="qheader,pident,evalue,bits,qcov,theader,tstart,tend"
    diamond_format="qseqid pident evalue bitscore qcovhsp stitle sstart send"
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

#----------------------#
# Calculate Alignments #
#----------------------#
touch aligner.log

if [ -s small.fasta ]; then
    if [ "$aliTool" == "blast" ]; then
        echo `blastp -version` >> aligner.log
        makeblastdb -in "$fasta" -dbtype prot -out orgdb >> aligner.log
        blastp -db orgdb -query small.fasta -qcov_hsp_perc $covcutoff -num_threads $n_threads -outfmt "6 $blast_format" > out.tsv
    fi

    if [ "$aliTool" == "diamond" ]; then
        echo `diamond --version` >> aligner.log
        diamond makedb --in "$fasta" -d orgdb >> aligner.log 2>&1
        diamond blastp -d orgdb.dmnd -q small.fasta  \
          --threads $n_threads \
          --out out.tsv \
          --outfmt 6 $diamond_format \
          --query-cover $covcutoff >> aligner.log 2>&1
    fi

    if [ "$aliTool" == "mmseqs2" ]; then
        echo `mmseqs version` >> aligner.log
        mmseqs createdb "$fasta" targetDB >> aligner.log 2>&1
        mmseqs createdb small.fasta queryDB >> aligner.log 2>&1
        mmseqs search queryDB targetDB resultDB $tmpdir --threads $n_threads -c 0.$covcutoff >> aligner.log 2>&1
        mmseqs convertalis queryDB targetDB resultDB out.tsv \
          --format-output "$mmseqs_format" >> aligner.log 2>&1
          
        sed -Ei 's/^([^ ]+) [^\t]+/\1/' out.tsv # get the fastq sequence identifier from query header (everything between the leading ">" and the first space).
    fi
else
    touch out.tsv
fi

#----------------------#
# Analyse Alignments   #
#----------------------#

Rscript $dir/analyse_alignments_transport.R $bitcutoff $identcutoff $nouse_alternatives $verbose


#-------------------------#
# Exporting results table # 
#-------------------------#

# add gapseq vesion and sequence database status to table comments head
gapseq_version=$($dir/.././gapseq -v | head -n 1)
seqdb_version=$(md5sum $dir/../dat/seq/transporter.fasta | cut -c1-7)
seqdb_date=$(stat -c %y $dir/../dat/seq/transporter.fasta | cut -c1-10)

if [ $input_mode == "nucl" ] && [ $newtranslate == "true" ]; then
    gzip -c $fasta > "$output_dir/${fastaID}.faa.gz"
    mv ${fastaID}.gff "$output_dir/${fastaID}.gff"
fi

genome_info="genome_format=${input_mode}"
[[ $input_mode == "nucl" ]] && faamd5=`md5sum $output_dir/${fastaID}.faa.gz | cut -c1-7` && genome_info="${genome_info};translation_md5=${faamd5}"
[[ $input_mode == "nucl" ]] && [[ $newtranslate == "true" ]] && genome_info="${genome_info};translation_table=${transl_table}"

sed -i "1s/^/# $gapseq_version\n/" transporter.tbl
sed -i "2s/^/# Transporter sequence DB md5sum: $seqdb_version ($seqdb_date)\n/" transporter.tbl
sed -i "3s/^/# $genome_info\n/" transporter.tbl

cp transporter.tbl $output_dir/${fastaID}-Transporter.tbl
cp aligner.log $output_dir/${fastaID}-Transporter-find_aligner.log

# finishing
end_time=`date +%s`
echo Running time: `expr $end_time - $start_time` s.

