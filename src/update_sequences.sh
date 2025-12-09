#!/bin/bash
zenodoID=10047603
zenodoRecord=16908828
taxonomy=Bacteria

quite=false

# paths
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
export LC_NUMERIC="en_US.UTF-8"
seqdb=$dir/../dat/seq
seqpath=$seqdb/$taxonomy
userdir=false
update_check=false
dbversion=$zenodoRecord

usage()
{
    echo "gapseq - update reference sequence database (update-sequences)"
    echo ""
    echo "Usage"
    echo "$0"
    echo "  -t Taxonomy of reference sequence to be updated. (default: $taxonomy)"
    echo "  -D Directory were the database will be stored/updated. (default: $seqdb)"
    echo "  -Z Database version as Zenodo ID (default: $zenodoRecord)"
    echo "  -c Check if a more recent sequence database is available on Zenodo. If option '-c' is used, option '-Z' is ignored and has no effect."
    echo "  -q No message output to prompt (quiet)."
    echo ""
    echo "Details:"
    echo "\"-Z\": This option expects the word 'latest' for the latest database version or the Zenodo record ID. All database records can be found here: https://doi.org/10.5281/zenodo.10047603 . The record ID is the last number in the DOI following the pattern 'zenodo.'"
exit 1
}

OPTIND=1         # Reset in case getopts has been used previously in the shell.
while getopts "h?t:D:Z:cq" opt; do
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
    c)
        update_check=true
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
    if [[ ! -w "$seqdb" ]] && [[ "$update_check" == false ]]; then
        echo "Error: No write permission for directory '$seqdb'. Contact your system administrator or choose a dabase directory, where you have write permission." >&2
        exit 1
    fi
    [[ $quite == false ]] && echo "Using custom datbase directory: $seqdb"
else
    # no -D provided → check default (<gapseq_dir>/dat/seq)
    if [[ ! -w "$seqdb" ]]; then
        # try fallback ~/.gapseq/seq
        if [[ -f "$seqdb/$taxonomy/version_seqDB.json" ]]; then
            echo "Warning: The default sequence database path '$seqdb' contains sequences for $taxonomy, but the directory is not writable."
            echo "         Will use the fallback directory '$HOME/.gapseq/seq' to store the new database."
            echo "         Remember to use the argument '-D $HOME/.gapseq/seq' when running 'gapseq find'."
        fi
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

# --- New Seq-DB available ? ---
if [[ "$update_check" == true ]]; then
    # get version info from current database
    current_version="NA"
    if [ -f "$seqpath/version_seqDB.json" ]; then
        current_version=`Rscript $dir/parse_seqDBjson.R "$seqpath/version_seqDB.json" zenodoID`
        current_versionNum=`Rscript $dir/parse_seqDBjson.R "$seqpath/version_seqDB.json" version`
        current_versionDate=`Rscript $dir/parse_seqDBjson.R "$seqpath/version_seqDB.json" date`
    fi

    if [[ "$current_version" == "NA" ]]; then
        echo "No current $taxonomy sequence database version found."
        echo "To download the gapseq version-specific sequence database, run:"
        echo "    gapseq update-sequences -t $taxonomy -D $seqdb"
        echo "To download the latest sequence database, run:"
        echo "    gapseq update-sequences -t $taxonomy -D $seqdb -Z latest"
        exit 0
    fi

    echo "Local $taxonomy sequence database: $current_versionNum (zenodoID: $current_version, date: $current_versionDate)"

    # check if newer database is on zenodo
    url_zenodorecord=https://zenodo.org/api/records/$zenodoID
    url_zenodocurrent=$(curl -Ls -o /dev/null -w %{url_effective} $url_zenodorecord)
    latest_version=$(basename "$url_zenodocurrent")

    if [[ "$latest_version" !=  "$current_version" ]]; then
        zen_record_summary=$(mktemp)
        wget -q -O "$zen_record_summary" $url_zenodocurrent
        zencommunity=$(Rscript -e 'args <- commandArgs(trailingOnly = TRUE); zenjs <- jsonlite::read_json(args[1]); cat(zenjs$metadata$communities[[1]]$id)' "$zen_record_summary")
        [[ "$zencommunity" != "gapseq" ]] && echo "Provided Zenodo-ID is not a valid gapseq record." && exit 1
        Rscript -e 'args <- commandArgs(trailingOnly = TRUE); zenjs <- jsonlite::read_json(args[1]); cat(jsonlite::toJSON(list(zenodoID = zenjs$id, version = zenjs$metadata$version, date = zenjs$created)))' "$zen_record_summary" > ${zen_record_summary}_2
        latest_version=`Rscript $dir/parse_seqDBjson.R "${zen_record_summary}_2" zenodoID`
        latest_versionNum=`Rscript $dir/parse_seqDBjson.R "${zen_record_summary}_2" version`
        latest_versionDate=`Rscript $dir/parse_seqDBjson.R "${zen_record_summary}_2" date`
        rm $zen_record_summary
        rm ${zen_record_summary}_2

        if [[ "$latest_version" == "NA" ]]; then
            "Error: Failed to fetch info for latest sequence database on zenodo."
            exit 1
        fi

        echo "[NOTE] A newer $taxonomy sequence database version exists: $latest_versionNum (zenodoID: $latest_version, date: $latest_versionDate)"
        echo "[NOTE] To update your local database run:"
        echo "[NOTE]     gapseq update-sequences -t $taxonomy -D $seqdb -Z latest"
    else
        echo "Reference sequences are up-to-date."
    fi
    exit 0
fi

if [[ "$dbversion" != "latest" ]]; then
  zenodoID=$dbversion
fi

[[ $quite == false ]] && echo "Checking sequence database for $taxonomy (DB-Path: $seqdb)"

url_zenodorecord=https://zenodo.org/api/records/$zenodoID
url_zenodocurrent=$(curl -Ls -o /dev/null -w %{url_effective} $url_zenodorecord)
url_zenodoseqs=$url_zenodocurrent/files/md5sums.txt/content
status_zenodo=$(curl -s --head -w %{http_code} $url_zenodoseqs -o /dev/null)

if [[ $status_zenodo -ge 500 ]]; then
    echo "Error: Failed to fetch info for latest sequence database on zenodo."
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

# get zenodo record info
zen_record_summary=$(mktemp)
wget -q -O "$zen_record_summary" $url_zenodocurrent
zencommunity=$(Rscript -e 'args <- commandArgs(trailingOnly = TRUE); zenjs <- jsonlite::read_json(args[1]); cat(zenjs$metadata$communities[[1]]$id)' "$zen_record_summary")
[[ "$zencommunity" != "gapseq" ]] && echo "Provided Zenodo-ID is not a valid gapseq record." && exit 1
Rscript -e 'args <- commandArgs(trailingOnly = TRUE); zenjs <- jsonlite::read_json(args[1]); cat(jsonlite::toJSON(list(zenodoID = zenjs$id, version = zenjs$metadata$version, date = zenjs$created)))' "$zen_record_summary" > ${zen_record_summary}_2
zen_version=`Rscript $dir/parse_seqDBjson.R "${zen_record_summary}_2" zenodoID`
zen_versionNum=`Rscript $dir/parse_seqDBjson.R "${zen_record_summary}_2" version`
zen_versionDate=`Rscript $dir/parse_seqDBjson.R "${zen_record_summary}_2" date`
rm $zen_record_summary
rm ${zen_record_summary}_2

if [[ $update_req -eq 0 ]]; then
  echo "Reference sequences for $taxonomy are already the desired version: $zen_versionNum (zenodoID: $zen_version, date: $zen_versionDate)"
  exit 0
fi

# Download and extract new sequences
if [[ $update_req -eq 1 ]]; then
  echo "Updating $taxonomy reference sequences to version: $zen_versionNum (zenodoID: $zen_version, date: $zen_versionDate)"

  # download+extract taxonomy-specific archive
  allseqs=$(mktemp)
  wget -q -O "$seqdb/$taxonomy.tar.gz" $url_zenodocurrent/files/$taxonomy.tar.gz/content
  tar xzf $seqdb/$taxonomy.tar.gz -C $seqdb/
  rm $seqdb/$taxonomy.tar.gz

  # extract rev/unrev/rxn archives
  tar xzf $seqpath/rev/sequences.tar.gz -C $seqpath/rev/
  tar xzf $seqpath/unrev/sequences.tar.gz -C $seqpath/unrev/
  tar xzf $seqpath/rxn/sequences.tar.gz -C $seqpath/rxn/

  # create/update zenodo record stamp
  zen_record_summary=$(mktemp)
  wget -q -O "$zen_record_summary" $url_zenodocurrent
  Rscript -e 'args <- commandArgs(trailingOnly = TRUE); zenjs <- jsonlite::read_json(args[1]); cat(jsonlite::toJSON(list(zenodoID = zenjs$id, version = zenjs$metadata$version, date = zenjs$created)))' "$zen_record_summary" > $seqpath/version_seqDB.json
  rm $zen_record_summary

  echo "Reference sequences updated."
fi

rm $zen_sums.$taxonomy

