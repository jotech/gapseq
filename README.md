

# gapseq
_Informed prediction and analysis of bacteria metabolic pathways and genome-scale networks_

_gapseq_ is designed to combine metabolic pathway analysis with metabolic network reconstruction and curation.
Based on genomic information and databases for pathways and reactions, _gapseq_ can be used for:
- prediction of metabolic pathways from various databases
- transporter inference
- metabolic model creation
- multi-step gap filling 

_Cite gapseq_ (preprint):\
gapseq: Informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models.\
Johannes Zimmermann, Christoph Kaleta, Silvio Waschina.\
_bioRxiv_ 2020.03.20.000737; doi: https://doi.org/10.1101/2020.03.20.000737 


## Installation
### Ubuntu/Debian/Mint
```
sudo apt install ncbi-blast+ git libglpk-dev r-base-core exonerate bedtools barrnap
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

### Centos/Fedora/RHEL
```
sudo yum install ncbi-blast+ git glpk-devel BEDTools exonerate hmmer
git clone https://github.com/tseemann/barrnap.git
export PATH=$PATH:barrnap/bin/barrnap # needs to be permanent => .bashrc ?
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

### MacOS
using [homebrew](https://brew.sh)
```
brew install git glpk blast bedtools r brewsci/bio/barrnap
R -e 'install.packages(c("data.table", "stringr", "sybil", "getopt", "reshape2", "doParallel", "foreach", "R.utils", "stringi", "glpkAPI"))'
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")'
git clone https://github.com/jotech/gapseq && cd gapseq
```

### Troubleshooting
- at least ncbi **blast** version 2.2.27 (4/2013) is needed. If your disribution only contains an older version, try to download a binary directly from [ncbi](https://shorturl.at/jkAH0)
- If you are getting the installation error ``'lib = "../R/library"' is not writable`` while installing the **R packages**, then try this command beforehand:
```
Rscript -e 'if( file.access(Sys.getenv("R_LIBS_USER"), mode=2) == -1 ) dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)'
```
- glpkAPI install error ``checking for library containing glp_create_prob... no`` although libglpk has been installed, try:
```
wget https://cran.r-project.org/src/contrib/glpkAPI_1.3.2.tar.gz
R CMD INSTALL --configure-args="--enable-gmp=no" glpkAPI_1.3.2.tar.gz
```

## Quickstart
 This predicts network candidate reactions, builds a draft model and performs gap filling:
```
./gapseq doall toy/myb71.fna
```
Do the same but with a defined medium for gap filling:
```
./gapseq doall toy/ecoli.fna.gz dat/media/MM_glu.csv
```

## Pathway analysis
Search for chitin pathways:
```
./gapseq find -p chitin toy/myb71.fna.gz
```
Check for a certain enzyme availability, for example cytochrome c oxidases (oxidase test)
```
./gapseq find -e 1.9.3.1 toy/ecoli.fna.gz
```
Search for enzymes by name:
```
./gapseq find -r ligninase toy/myb71.fna.gz
```

## Creation and gap filling of metabolic models
1) Metabolic pathway analysis
```
./gapseq find -p all toy/myb71.fna
./gapseq find-transport toy/myb71.fna
```

2) Creation of draft model
```
 ./gapseq draft -r toy/myb71-all-Reactions.tbl -t toy/myb71-Transporter.tbl -p toy/myb71-all-Pathways.tbl -c toy/myb71.fna.gz
 ```

3) Gap filling
```
./gapseq fill -m toy/myb71-draft.RDS -c toy/myb71-rxnWeights.RDS -g toy/myb71-rxnXgenes.RDS -n dat/media/TSBmed.csv
```


------

###### LICENSE

Copyright 2020 Johannes Zimmermann, Christoph Kaleta, & Silvio Waschina; University of Kiel, Germany

GNU General Public License version 3.0 ([GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)) is applied to all copyrightable parts of gapseq. **gapseq** uses information on biochemical reactions, compounds, compartments, enzymes, and biological sequences from different external sources. The copyright and licensing terms for each of the resources are listed and cross-linked below. Identifiers for reactions, enzymes, compounds, and compartments may be identical to the external sources but can also differ to those. In both cases, the data from **gapseq** may be considered to be subject to the original copyright and licensing restrictions of the external resource.

- **MNXref**: Copyright 2011-2019 SystemsX, SIB Swiss Institute of Bioinformatics. 

  Licensed under a Creative Commons Attribution 4.0 International License.

  Link to license: https://creativecommons.org/licenses/by/4.0/

  Website: https://www.metanetx.org/

- **MetaCyc**: Copyright © SRI International 1999-2020, Marine Biological Laboratory  1998-2001, DoubleTwist Inc 1998-1999.  

  Link to license: https://metacyc.org/ptools-academic-license.shtml .

  Website: https://metacyc.org/

- **MODELSEED**: Copyright 2015 ModelSEED

  Licensed under Creative Commons  Attribution 4.0 International License.

  Link to license: https://creativecommons.org/licenses/by/4.0/ 

  Website: https://metacyc.org/

- **KEGG**: [Copyright](https://www.kegg.jp/kegg/legal.html) 1995-2020 [Kanehisa Laboratories](https://www.kanehisa.jp/)

  For license terms see file `dat/licenses/LICENSE.kegg` .

  Website: http://www.kegg.jp

- **BRENDA**: Copyright 2020 Prof. Dr. D. Schomburg, Technische Universität Braunschweig,  BRICS, Department of Bioinformatics and Biochemistry, Rebenring 56, 38106 Braunschweig, Germany.

  Licensed under the Creative Commons Attribution License [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) is applied to all[ copyrightable parts](https://wiki.creativecommons.org/wiki/Data#Can_databases_be_released_under_CC_licenses.3F) of BRENDA.

  Link to license: https://creativecommons.org/licenses/by/4.0/

  Website: https://www.brenda-enzymes.org/

- **UNIPROT**: Copyright 2002 –2020 [UniProt Consortium](https://www.uniprot.org/help/about)

  Licensed under the Creative Commons Attribution License [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) is applied to all[ copyrightable parts](https://wiki.creativecommons.org/wiki/Data#Can_databases_be_released_under_CC_licenses.3F) of UNIPROT.

  Link to license: https://creativecommons.org/licenses/by/4.0/

  Website: https://www.uniprot.org/

- **TCDB**: Copyright 2005 - 2020 Saier Lab.

  The text of the TCDB website (TCDB.ORG) is available for modification and reuse under  the terms of the Creative Commons Attribution-Sharealike 3.0 Unported  License and the GNU Free Documentation License.

  Link to license: https://creativecommons.org/licenses/by/4.0/

  Website: http://www.tcdb.org/

- **BIGG**: Copyright © 2019 The Regents of the University of California 

  All Rights Reserved by the licenser. For license terms and conditions see file `dat/licenses/LICENSE.bigg` .
  
  Website: http://bigg.ucsd.edu/
