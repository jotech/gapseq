#!/bin/bash
zenodoID=10047603
taxonomy=Bacteria

quite=false

# paths
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
export LC_NUMERIC="en_US.UTF-8"
seqdb=$dir/../dat/seq
seqpath=$seqdb/$taxonomy
dbversion=latest
userdir=false

usage()
{
    echo "gapseq - update reference sequence database (update-sequences)"
    echo ""
    echo "Usage"
    echo "$0"
    echo "  -t Taxonomy of reference sequence to be updated. (default: $taxonomy)"
    echo "  -D Directory were the database will be stored/updated. (default: $seqdb)"
    echo "  -Z Database version as Zenodo ID (default: $dbversion)"
    echo "  -q No message output to prompt (quiet)."
    echo ""
    echo "Details:"
    echo "\"-Z\": This option expects the the word 'latest' for the latest database version or the Zenodo record ID. All database records can be found here: https://doi.org/10.5281/zenodo.10047603 . The record ID is the last number in the DOI following the pattern 'zenodo.'"
exit 1
}

OPTIND=1         # Reset in case getopts has been used previously in the shell.
while getopts "h?t:D:Z:q" opt; do
    case "$opt" in
    h|\?)
        usage
        exit 0
        ;;
    t)
        taxonomy=$OPTARG
        ;;
    D)
        seqdb=$OPTARG
        userdir=true
        ;;
    Z)
        dbversion=$OPTARG
        ;;
    q)
        quite=true
        ;;
    esac
done

# --- Directory checks ---
if [[ "$userdir" == true ]]; then
    # user provided -D
    seqdb=$(readlink -f "$seqdb")
    mkdir -p "$seqdb" || {
        echo "Error: could not create database directory '$seqdb'." >&2
        exit 1
    }
    if [[ ! -w "$seqdb" ]]; then
        echo "Error: directory '$seqdb' is not writable." >&2
        exit 1
    fi
    [[ $quite == false ]] && echo "Using custom datbase directory: $seqdb"
else
    # no -D provided → check default
    if [[ ! -w "$seqdb" ]]; then
        # try fallback ~/.gapseq/seq
        seqdb="$HOME/.gapseq/seq"
        mkdir -p "$seqdb" || {
            echo "Error: could not create fallback directory '$seqdb'." >&2
            exit 1
        }
        [[ $quite == false ]] && echo "Using fallback directory: $seqdb"
    fi
fi

seqpath=$seqdb/$taxonomy
mkdir -p $seqpath/rev $seqpath/unrev $seqpath/rxn


if [[ "$dbversion" != "latest" ]]; then
  zenodoID=$dbversion
fi

echo "Checking sequence database for $taxonomy (DB-Path: $seqdb; version: $dbversion)"

url_zenodorecord=https://zenodo.org/records/$zenodoID
url_zenodocurrent=$(curl -Ls -o /dev/null -w %{url_effective} $url_zenodorecord)
url_zenodoseqs=$url_zenodocurrent/files/md5sums.txt
status_zenodo=$(curl -s --head -w %{http_code} $url_zenodoseqs -o /dev/null)

if [[ $status_zenodo -ge 500 ]]; then
    echo "No sequence archive found, manual download needed."
    exit 1
fi

dir_rev=$seqpath/rev
dir_unrev=$seqpath/unrev
dir_rxn=$seqpath/rxn

# get md5 checksums from the online data version
zen_sums=$(mktemp)
wget -q -O $zen_sums $url_zenodoseqs
cat $zen_sums | grep -w $taxonomy > $zen_sums.$taxonomy
rm $zen_sums

# check for mismatches of local md5sums and online version
update_req=0
while read -r md5sum filepath; do
  if [ ! -f "$seqdb/$filepath" ]; then
    update_req=1
    break
  fi

  # Calculate the MD5 sum of local file
  calculated_md5=$(md5sum "$seqdb/$filepath" | cut -d ' ' -f 1)

  # Compare the calculated MD5 sum with the one from the file
  if [ "$md5sum" != "$calculated_md5" ]; then
    update_req=1
    break
  fi
done < $zen_sums.$taxonomy

if [[ $update_req -eq 0 ]]; then
  echo "Reference sequences are already the desired version '$dbversion'."
  exit 0
fi

# Download and extract new sequences
if [[ $update_req -eq 1 ]]; then
  echo "Updating reference sequences to version '$dbversion' ..."

  # download+extract taxonomy-specific archive
  allseqs=$(mktemp)
  wget -q -O "$seqdb/$taxonomy.tar.gz" $url_zenodocurrent/files/$taxonomy.tar.gz
  tar xzf $seqdb/$taxonomy.tar.gz -C $seqdb/
  rm $seqdb/$taxonomy.tar.gz

  # extract rev/unrev/rxn archives
  tar xzf $seqpath/rev/sequences.tar.gz -C $seqpath/rev/
  tar xzf $seqpath/unrev/sequences.tar.gz -C $seqpath/unrev/
  tar xzf $seqpath/rxn/sequences.tar.gz -C $seqpath/rxn/

  echo "Reference sequences updated."
fi

rm $zen_sums.$taxonomy

