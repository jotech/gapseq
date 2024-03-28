# gapseq
_Informed prediction and analysis of bacterial metabolic pathways and genome-scale networks_  
[![Documentation Status](https://readthedocs.org/projects/gapseq/badge/?version=latest)](https://gapseq.readthedocs.io/en/latest/?badge=latest)
[![DOI:10.1186/s13059-021-02295-1](https://zenodo.org/badge/DOI/10.1186/s13059-021-02295-1.svg)](https://doi.org/10.1186/s13059-021-02295-1)  


_gapseq_ is designed to combine metabolic pathway analysis with metabolic network reconstruction and curation.
Based on genomic information and databases for pathways and reactions, _gapseq_ can be used for:
- prediction of metabolic pathways from various databases
- transporter inference
- metabolic model construction
- multi-step gap filling 

<p align="center">
<img src="https://github.com/jotech/gapseq/raw/master/docs/gfx/flowchart.png" alt="alt text" title="Title" width="600">
</p>

## Publication
Zimmermann, J., Kaleta, C. & Waschina, S. [gapseq: informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models.](https://doi.org/10.1186/s13059-021-02295-1) *Genome Biology* 22, 81 (2021)


## Installation
The latest release can be downloaded [here](https://github.com/jotech/gapseq/releases).
Besides this, the current development version can be accessed via:
```
git clone https://github.com/jotech/gapseq
```
Detailed information on [installation and troubleshooting](https://github.com/jotech/gapseq/blob/master/docs/install.md).


## Quickstart
For detailed use cases and full tutorials, see the [documentation](https://gapseq.readthedocs.io/).

Prediction of network candidate reactions, building of a draft model and gap filling:
```
./gapseq doall toy/myb71.fna
```
Do the same but with a defined medium for gap filling:
```
./gapseq doall toy/ecoli.fna.gz dat/media/MM_glu.csv
```


------

## LICENSE

Copyright 2020 Johannes Zimmermann, Christoph Kaleta, & Silvio Waschina; University of Kiel, Germany

GNU General Public License version 3.0 ([GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)) is applied to all copyrightable parts of gapseq. **gapseq** uses information on biochemical reactions, compounds, compartments, enzymes, and biological sequences from different external sources. The copyright and licensing terms for each of the resources are listed and cross-linked below. Identifiers for reactions, enzymes, compounds, and compartments may be identical to the external sources but can also differ to those. In both cases, the data from **gapseq** may be considered to be subject to the original copyright and licensing restrictions of the external resource.

- **MNXref**: Copyright 2011-2019 SystemsX, SIB Swiss Institute of Bioinformatics. 
  Licensed under a Creative Commons Attribution 4.0 International License.
  Link to license: https://creativecommons.org/licenses/by/4.0/
  Website: https://www.metanetx.org/
- **MetaCyc**: Copyright © SRI International 1999-2020, Marine Biological Laboratory  1998-2001, DoubleTwist Inc 1998-1999.  
  Link to license: https://metacyc.org/ptools-academic-license.shtml .
  Website: https://metacyc.org/
- **MODELSEED**: Copyright 2015 ModelSEED.
  Licensed under Creative Commons  Attribution 4.0 International License.
  Link to license: https://creativecommons.org/licenses/by/4.0/ 
  Website: https://modelseed.org/
- **KEGG**: [Copyright](https://www.kegg.jp/kegg/legal.html) 1995-2020 [Kanehisa Laboratories](https://www.kanehisa.jp/).
  For license terms see file `dat/licenses/LICENSE.kegg` .
  Website: http://www.kegg.jp
- **BRENDA**: Copyright 2020 Prof. Dr. D. Schomburg, Technische Universität Braunschweig,  BRICS, Department of Bioinformatics and Biochemistry, Rebenring 56, 38106 Braunschweig, Germany.
  Licensed under the Creative Commons Attribution License [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) is applied to all[ copyrightable parts](https://wiki.creativecommons.org/wiki/Data#Can_databases_be_released_under_CC_licenses.3F) of BRENDA.
  Link to license: https://creativecommons.org/licenses/by/4.0/
  Website: https://www.brenda-enzymes.org/
- **UNIPROT**: Copyright 2002 –2020 [UniProt Consortium](https://www.uniprot.org/help/about).
  Licensed under the Creative Commons Attribution License [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) is applied to all[ copyrightable parts](https://wiki.creativecommons.org/wiki/Data#Can_databases_be_released_under_CC_licenses.3F) of UNIPROT.
  Link to license: https://creativecommons.org/licenses/by/4.0/
  Website: https://www.uniprot.org/
- **TCDB**: Copyright 2005 - 2020 Saier Lab.
  The text of the TCDB website (TCDB.ORG) is available for modification and reuse under  the terms of the Creative Commons Attribution-Sharealike 3.0 Unported  License and the GNU Free Documentation License.
  Link to license: https://creativecommons.org/licenses/by-sa/3.0/de/
  Website: http://www.tcdb.org/
- **BIGG**: Copyright © 2019 The Regents of the University of California. 
  All Rights Reserved by the licenser. For license terms and conditions see file `dat/licenses/LICENSE.bigg` .
  Website: http://bigg.ucsd.edu/
- **GTDBtk**: Copyright 2017 Pierre-Alain Chaumeil.
  Licensed under the [GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) license.
  Link to license: https://github.com/Ecogenomics/GTDBTk/blob/master/LICENSE
  Website: https://ecogenomics.github.io/GTDBTk/

## Citation
- Zimmermann, J., Kaleta, C. & Waschina, S. gapseq: informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models. Genome Biol 22, 81 (2021). https://doi.org/10.1186/s13059-021-02295-1
- De Bernardini, N., Zampieri, G., Campanaro, S., Zimmermann, J., Waschina, S. & Treu, L. pan-Draft: Automated reconstruction of species-representative metabolic models from multiple genomes. *Currently under revision*
