The biomass reaction is an essential part in Flux Balance Analysis of microbial growth. The composition of the biomass reaction (i.e. stoichiometry of  molecular cell constituents) should ideally reflect the molecular makeup of the organism of interest.

#### Biomass reaction presets

gapseq currently provides three template biomass reactions: One for archaea and two for bacteria (Gram positive and Gram negative). 

|                                          | Bacteria (Gram negative)                                     | Bacteria (Gram positive)                                     | Archaea                                                      |
| ---------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| gapseq parameter for draft construction* | `-b Gram_neg`                                                | `-b Gram_pos`                                                | `-b archaea`                                                 |
| Biomass reaction file                    | `biommass_Gram_neg.json`                                     | `biommass_Gram_pos.json`                                     | `biommass_archaea.json`                                      |
|                                          |                                                              |                                                              |                                                              |
| <u>*Originally derived from:*</u>        |                                                              |                                                              |                                                              |
| Reference organism                       | *Escherichia coli*                                           | *Bacillus subtilis*                                          | *Methanosarcina barkeri*                                     |
| Reference publication**                  | [Orth *et al.* (2011) *Mol Syst Biol*](https://doi.org/10.1038/msb.2011.65) | [Oh *et al.* (2007) *J Biol Chem*](https://doi.org/10.1074/jbc.m703759200) | [Feist *et al.* (2006) *Mol Syst Biol*](https://doi.org/10.1038/msb4100046) |
|                                          |                                                              |                                                              |                                                              |
| *<u>Constituents groups</u>* in % (w/w ) |                                                              |                                                              |                                                              |
| DNA                                      | 3.1 %                                                        | 2.60 %                                                       | 4 %                                                          |
| RNA                                      | 20.5 %                                                       | 6.55 %                                                       | 24  %                                                        |
| Protein                                  | 55 %                                                         | 52.84 %                                                      | 63 %                                                         |
| Inorganic compounds                      | 0.5 %                                                        | 0.50 %                                                       | -                                                            |
| Co-factors & soluble compounds           | 3.4 %                                                        | 4.45 %                                                       | 4 %                                                          |
| Lipids                                   | 9.1 %                                                        | 7.60 %                                                       | 4.9 %                                                        |
| Cell wall                                | 8.4 %                                                        | 25.46 %                                                      | -                                                            |
| Polysaccharides                          | -                                                            | -                                                            | 0.1 %                                                        |
|                                          |                                                              |                                                              |                                                              |

\* If no biomass reaction is defined using this command, gapseq will predict the most likely biomass reaction based on the 16S rRNA gene sequence using a pre-trained classifier.

\** When using gapseq for publications, please also cite these original research papers to acknowledge their work. And of course gapseq :)



#### User-defined biomass reactions

We are currently working on an additional feature to gapseq that allows users to specify biomass reactions that are different from the above mentioned presets. We will describe how to use this feature here as soon as it is implemented.
