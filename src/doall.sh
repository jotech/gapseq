#!/bin/bash

bitcutoff=200
min_bs_core=100
identcutoff=0   # cutoff blast: min identity
covcutoff=75 # cutoff blast: min coverage
output_dir=.
user_temp=false
aliTool="blast"
medium="auto"
taxonomy="auto"
user_temp_opt=""
verbose=0

OS=$(uname -s)
if [ "$OS" = "Darwin" -o "$OS" = "FreeBSD" ]; then
	n_threads=$(sysctl hw.ncpu|cut -f2 -d' ')
else
	n_threads=`grep -c ^processor /proc/cpuinfo`
fi

usage()
{
    echo "gapseq - perfoming all reconstruction steps (doall)"
    echo ""
    echo "Usage"
    echo "$0 file.fasta"
    echo "  -b bit score cutoff for local alignment (default: $bitcutoff)"
    echo "  -i identity cutoff for local alignment (default: $identcutoff)"
    echo "  -c coverage cutoff for local alignment (default: $covcutoff)"
    echo "  -l Reactions with an associated blast-hit with a bitscore below this value will be considered just as reactions that have no blast hit."
    echo "  -t Taxonomic range for reference sequences to be used. (Bacteria, Archaea, auto; default: $taxonomy)."
    echo "  -m tab- or komma separated table for media components. Requires three named columns: 1 - \"compounds\" (for metab. IDs), 2 - \"name\" (metab. name), 3 - \"maxFlux\" (maximum inflow flux). If \"auto\" (the default), a growth medium will be predicted."
    echo "  -f Path to directory, where output files will be saved (default: current directory)"
    echo "  -v Verbose level, 0 for nothing, 1 for full (default $verbose)"
    echo "  -T Set user-defined temporary folder (default: $user_temp)"
    echo "  -K Number of threads. If option is not provided, number of available CPUs will be automatically determined."
    echo "  -A Tool to be used for sequence alignments (blast, mmseqs2, diamond; default: $aliTool)"
    echo ""
exit 1
}

OPTIND=1         # Reset in case getopts has been used previously in the shell.
while getopts "h?b:i:c:l:t:m:f:v:T:K:A:" opt; do
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
    l)
        min_bs_core=$OPTARG
        ;;
    t)  
        taxonomy=$OPTARG
        ;;
    m)
        medium=$OPTARG
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
        user_temp_opt="-T ${user_temp_folder}"
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

file=$(readlink -f $1)
base=$(basename "$file")
id="${base%.*}"
[[ ! -s "$file" ]]  && usage
[[ $file == *.gz ]] && id="${id%.*}" 

$dir/gapseq_find.sh -v $verbose -b $bitcutoff -p all -t $taxonomy -K $n_threads -A $aliTool -f $output_dir $user_temp_opt "$file"
$dir/transporter.sh -v $verbose -b $bitcutoff -K $n_threads -A $aliTool -f $output_dir $user_temp_opt "$file"
Rscript $dir/generate_GSdraft.R -r "$output_dir/$id-all-Reactions.tbl" -t "$output_dir/$id-Transporter.tbl" -u $bitcutoff -l $min_bs_core -p "$output_dir/$id-all-Pathways.tbl" -b $taxonomy -f $output_dir
if [ $medium == "auto" ]; then
    Rscript $dir/predict_medium.R -m "$output_dir/${id}-draft.RDS" -p "$output_dir/$id-all-Pathways.tbl" -f $output_dir
    medium="$output_dir/${id}-medium.csv"
fi
Rscript $dir/gf.suite.R -m "$output_dir/${id}-draft.RDS" -n "$medium" -c "$output_dir/${id}-rxnWeights.RDS" -b $min_bs_core -g "$output_dir/${id}-rxnXgenes.RDS" -f $output_dir


