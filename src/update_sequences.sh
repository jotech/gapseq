#!/bin/bash

# taxonomy
if [[ -z "$1" ]]; then
    taxonomy="Bacteria"
else taxonomy=$1
fi

# paths
curdir=$(pwd)
path=$(readlink -f "$0")
dir=$(dirname "$path")
export LC_NUMERIC="en_US.UTF-8"
seqpath=$dir/../dat/seq/$taxonomy
seqpath_user=$dir/../dat/seq/$taxonomy/user
mkdir -p $seqpath/rev $seqpath/unrev $seqpath_user $seqpath/rxn

echo Checking updates for $taxonomy $seqpath

url_rev=ftp://ftp.rz.uni-kiel.de/pub/medsystbio/$taxonomy/rev/sequences.tar.gz
url_unrev=ftp://ftp.rz.uni-kiel.de/pub/medsystbio/$taxonomy/unrev/sequences.tar.gz
url_rxn=ftp://ftp.rz.uni-kiel.de/pub/medsystbio/$taxonomy/rxn/sequences.tar.gz

status_rev=$(curl -s --head -w %{http_code} $url_rev -o /dev/null)
status_unrev=$(curl -s --head -w %{http_code} $url_unrev -o /dev/null)
status_rxn=$(curl -s --head -w %{http_code} $url_rxn -o /dev/null)

if [[ $status_rev -eq 550 ]] && [[ $status_unrev -eq 550 ]] && [[ $status_rxn -eq 550 ]]; then
    echo "No sequence archive found, manual download needed."
    exit 1
fi

dir_rev=$seqpath/rev
dir_unrev=$seqpath/unrev
dir_rxn=$seqpath/rxn

# check modification time
if [[ -s $dir_rev/sequences.tar.gz ]]; then
    mod_rev=$(stat -c %Y $dir_rev/sequences.tar.gz)
else
    mod_rev=0
fi
if [[ -s $dir_unrev/sequences.tar.gz ]]; then
    mod_unrev=$(stat -c %Y $dir_unrev/sequences.tar.gz)
else
    mod_unrev=0
fi
if [[ -s $dir_rxn/sequences.tar.gz ]]; then
    mod_rxn=$(stat -c %Y $dir_rxn/sequences.tar.gz)
else
    mod_rxn=0
fi

# download if newer
cd $dir_rev && wget -nv -N $url_rev
cd $dir_unrev && wget -nv -N $url_unrev
cd $dir_rxn && wget -nv -N $url_rxn

# extract if new files were downloaded
modnew_rev=$(stat -c %Y $dir_rev/sequences.tar.gz)
modnew_unrev=$(stat -c %Y $dir_unrev/sequences.tar.gz)
modnew_rxn=$(stat -c %Y $dir_rxn/sequences.tar.gz)

if [[ $modnew_rev -gt $mod_rev ]]; then
    tar xzf $seqpath/rev/sequences.tar.gz -C $seqpath/rev/
    echo "Updated $taxonomy reviewed sequences"
else
    echo "$taxonomy reviewed sequences already up-to-date"
fi
if [[ $modnew_unrev -gt $mod_unrev ]]; then
    tar xzf $seqpath/unrev/sequences.tar.gz -C $seqpath/unrev/
    echo "Updated $taxonomy unreviewed sequences"
else
    echo "$taxonomy unreviewed sequences already up-to-date"
fi
if [[ $modnew_rxn -gt $mod_rxn ]]; then
    tar xzf $seqpath/rxn/sequences.tar.gz -C $seqpath/rxn/
    echo "Updated $taxonomy additional reaction sequences"
else
    echo "$taxonomy additional reaction sequences already up-to-date"
fi
