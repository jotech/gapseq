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
mkdir -p $seqpath/rev $seqpath/unrev $seqpath_user

echo Checking updates for $taxonomy $seqpath

url_bacrev=ftp://ftp.rz.uni-kiel.de/pub/medsystbio/Bacteria/rev/sequences.tar.gz
url_bacunrev=ftp://ftp.rz.uni-kiel.de/pub/medsystbio/Bacteria/unrev/sequences.tar.gz

dir_bacrev=$seqpath/rev
dir_bacunrev=$seqpath/unrev

# check modification time
if [[ -s $dir_bacrev/sequences.tar.gz ]]; then
    mod_bacrev=$(stat -c %Y $dir_bacrev/sequences.tar.gz)
else
    mod_bacrev=0
fi
if [[ -s $dir_bacunrev/sequences.tar.gz ]]; then
    mod_bacunrev=$(stat -c %Y $dir_bacunrev/sequences.tar.gz)
else
    mod_bacunrev=0
fi

# download if newer
cd $dir_bacrev && wget -nv -N $url_bacrev
cd $dir_bacunrev && wget -nv -N $url_bacunrev

# extract if new files were downloaded
modnew_bacrev=$(stat -c %Y $dir_bacrev/sequences.tar.gz)
modnew_bacunrev=$(stat -c %Y $dir_bacunrev/sequences.tar.gz)

if [[ $modnew_bacrev -gt $mod_bacrev ]]; then
    tar xzf $seqpath/rev/sequences.tar.gz -C $seqpath/rev/
    echo "Updated bacteria reviewed sequences"
else
    echo "Bacteria reviewed sequences already up-to-date"
fi
if [[ $modnew_bacunrev -gt $mod_bacunrev ]]; then
    tar xzf $seqpath/unrev/sequences.tar.gz -C $seqpath/unrev/
    echo "Updated bacteria unreviewed sequences"
else
    echo "Bacteria unreviewed sequences already up-to-date"
fi
