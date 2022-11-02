# Reference HMM protein family models

## Prediction of domain (Bacteria/Archaea) via HMM profiles

The files stored in this directory contains the required data to infer the organism's domain (Bacteria/Archaea) from the input genome sequence. To this end, the file `domain.hmm.gz` contains HMM profiles of orthologous protein-coding genes (amino acid sequence) that are either indicative for Bacteria or Archaea. The assignment of HMM-profile to the specific domain is retrieved from [GTDBtk](https://github.com/Ecogenomics/GTDBTk) (package version R207_v2). 

- **GTDBtk**: Copyright 2017 Pierre-Alain Chaumeil. LICENSE: GPL-3.0. Resource website: https://ecogenomics.github.io/GTDBTk/

## Prediction of gram staining (pos/neg) via HMM profiles

Gram staining-associated protein HMM profiles were retrieved from the ["Protein Family Models" database](https://www.ncbi.nlm.nih.gov/protfam) from NCBI on November 2nd 2022. To this end, following query was used to find <u>Gram-positive</u> indicating profiles:

`(Gram-positive) NOT Gram-negative NOT response to Gram-positive`

Accordingly, <u>Gram-negative</u> indicating profiles via:

`(Gram-negative) NOT Gram-positive NOT response to Gram-negative`

The IDs of found HMM profiles and their Gram-staining association is provided in the file `gram_markers.tsv` and the HMM profiles in `gram.hmm.gz`.
