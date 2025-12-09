# Reference sequences

gapseq relies on reference amino acid sequences of proteins with an annotated enzymatic function. These sequences are obtained from UniProt.

We provide regularly updated [reference sequence databases for gapseq via Zenodo](https://doi.org/10.5281/zenodo.10047603). These databases can be installed and updated by:

```sh
gapseq update-sequences -t Bacteria # for bacteria
gapseq update-sequences -t Archaea # for archaea
```

By default, the references sequences are stored in `<gapseq_dir>/dat/seq/`, where "<gapseq_dir>" is the gapseq installation directory.

Users can also specify an alternative location in which to store the reference sequence database. This is particularly important when the user does not have write permissions in the gapseq installation directory.

To specify a non-default directory for the storage of the sequence database, use the option `-D` in `gapseq find` and `gapseq update-sequences`. E.g.:

```sh
# Download the database to a user-directory
gapseq update-sequences -t Bacteria -D $HOME/gapseq_seqDB

# Use the sequence database from the user directory in 'find'
gapseq find -p all -D $HOME/gapseq_seqDB -A diamond ecoli.faa.gz
```

### Fallback location for reference sequence database

If the gapseq installation directory is not writable, i.e. if the user does not have write permission, gapseq tries to save the sequence database in ~/.gapseq/seq. Users should therefore add the argument `-D ~/.gapseq/seq` when using the `gapseq find` module.

If the gapseq installation directory is not writable but contains a sequence database (for example, if it was installed by an administrator), the `gapseq find` module will use this database if no other path is specified using the `-D` option.

### 'gapseq version-specific reference sequence database' vs. 'latest reference sequence database'

From version 2.0 onwards, each version of gapseq (and each commit) is linked to a specific version of the sequence database. Therefore, the command `gapseq update-sequences -t Bacteria` may not retrieve the latest available sequence database on Zenodo. This behaviour ensures reproducibility of results from a gapseq version with the corresponding sequence database with which the gapseq version was tested.

To update the sequence DB to the latest version, simply run:

```sh
gapseq update-sequences -t Bacteria -Z latest`
```

Alternatively, update your gapseq version to the latest development version, which is usually linked to the latest sequence DB.

To install a specific database version run:

```sh
gapseq update-sequences -t Bacteria -Z <zenodoID>`
```

where `<zenodoID>` is a Zenodo record ID of the desired sequence database version. A list of all sequence database releases can be found [here](https://doi.org/10.5281/zenodo.10047603).

