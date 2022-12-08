# Biomass Reactions

The biomass reaction is an essential part in Flux Balance Analysis of microbial growth. The composition of the biomass reaction (i.e. stoichiometry of molecular cell constituents) should ideally reflect the molecular makeup of the organism of interest.

## Biomass reaction presets

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

\* If no biomass reaction is defined using this command, gapseq will predict the most likely biomass reaction based on the 16S rRNA gene sequence using a pre-trained classifier (when input is nucleotide fasta) or specific protein-coding marker genes (when input is protein fasta).

\** When using gapseq for publications, please also cite these original research papers to acknowledge their work. And of course gapseq :)

## User-defined biomass reactions

As an advanced usage feature, gapseq allows user-defined biomass reactions.

To use this feature, simple supply the path to your biomass definition JSON-file with the option `-b` for the `gapseq draft` module. For example:

```sh
gapseq draft -r genome-all-Reactions.tbl -t genome-Transporter.tbl -b ~/path/to/biomass/user_biomass.json -c genome.fna.gz -p genome-all-Pathways.tbl
```

The JSON-file for the biomass formulation needs to follow a specific format. An example for a user-defined biomass json-file:

```json
{ "id" : "user_biomass",
  "name" : "A user biomass definition",
  "ref"  : "Motivated by https://github.com/jotech/gapseq/issues/156",
  "energy_GAM" :  40,
  "domain" : "Bacteria",
  "met_groups" : [
    {
      "group_name" : "DNA",
      "mass" : 0.03,
      "unit_group" : "g",
      "unit_components" : "MOLFRACTION",
      "components" : [
        
        {
          "id"   : "cpd00115",
          "name" : "dATP",
          "comp" : "c",
          "coef" : 0.25,
          "link" : "cpd00012:-1"
        },
        
        {
          "id"   : "cpd00357",
          "name" : "dTTP",
          "comp" : "c",
          "coef" : 0.25,
          "link" : "cpd00012:-1"
        },
        
        {
          "id"   : "cpd00241",
          "name" : "dGTP",
          "comp" : "c",
          "coef" : 0.25,
          "link" : "cpd00012:-1"
        },
        
        {
          "id"   : "cpd00356",
          "name" : "dCTP",
          "comp" : "c",
          "coef" : 0.25,
          "link" : "cpd00012:-1"
        }
        
      ]
    },
    
    {
      "group_name" : "RNA",
      "mass" : 0.2,
      "unit_group" : "g",
      "unit_components" : "MOLFRACTION",
      "components" : [
        
        {
          "id"   : "cpd00002",
          "name" : "ATP",
          "comp" : "c",
          "coef" : 0.25,
          "link" : "cpd00012:-1"
        },
        
        {
          "id"   : "cpd00062",
          "name" : "UTP",
          "comp" : "c",
          "coef" : 0.25,
          "link" : "cpd00012:-1"
        },
        
        {
          "id"   : "cpd00038",
          "name" : "GTP",
          "comp" : "c",
          "coef" : 0.25,
          "link" : "cpd00012:-1"
        },
        
        {
          "id"   : "cpd00052",
          "name" : "CTP",
          "comp" : "c",
          "coef" : 0.25,
          "link" : "cpd00012:-1"
        }
        
      ]
    },
    
    {
      "group_name" : "Protein",
      "mass" : 0.55,
      "unit_group" : "g",
      "unit_components" : "MOLFRACTION",
      "components" : [
        
        {
          "id"   : "cpd00035",
          "name" : "L-alanine",
          "comp" : "c",
          "coef" : 0.1,
          "link" : "cpd00001:-1"
        },
        
        {
          "id"   : "cpd00051",
          "name" : "L-arginine",
          "comp" : "c",
          "coef" : 0.05,
          "link" : "cpd00001:-1"
        },
        
        {
          "id"   : "cpd00132",
          "name" : "L-asparagine",
          "comp" : "c",
          "coef" : 0.05,
          "link" : "cpd00001:-1"
        },
        
        {
          "id"   : "cpd00041",
          "name" : "L-Aspartate",
          "comp" : "c",
          "coef" : 0.05,
          "link" : "cpd00001:-1"
        },        
        
        {
          "id"   : "cpd00084",
          "name" : "L-Cysteine",
          "comp" : "c",
          "coef" : 0.01,
          "link" : "cpd00001:-1"
        },
        
        {
          "id"   : "cpd00053",
          "name" : "L-Glutamine",
          "comp" : "c",
          "coef" : 0.05,
          "link" : "cpd00001:-1"
        },        
        
        {
          "id"   : "cpd00023",
          "name" : "L-Glutamate",
          "comp" : "c",
          "coef" : 0.06,
          "link" : "cpd00001:-1"
        },        
        
        {
          "id"   : "cpd00033",
          "name" : "Glycine",
          "comp" : "c",
          "coef" : 0.1,
          "link" : "cpd00001:-1"
        },        
        
        {
          "id"   : "cpd00119",
          "name" : "L-Histidine",
          "comp" : "c",
          "coef" : 0.02,
          "link" : "cpd00001:-1"
        },
        
        {
          "id"   : "cpd00322",
          "name" : "L-Isoleucine",
          "comp" : "c",
          "coef" : 0.05,
          "link" : "cpd00001:-1"
        },        
        
        {
          "id"   : "cpd00107",
          "name" : "L-Leucine",
          "comp" : "c",
          "coef" : 0.09,
          "link" : "cpd00001:-1"
        },        
        
        {
          "id"   : "cpd00039",
          "name" : "L-Lysine",
          "comp" : "c",
          "coef" : 0.06,
          "link" : "cpd00001:-1"
        },
        
        {
          "id"   : "cpd00060",
          "name" : "L-Methionine",
          "comp" : "c",
          "coef" : 0.02,
          "link" : "cpd00001:-1"
        },        
        
        {
          "id"   : "cpd00066",
          "name" : "L-Phenylalanine",
          "comp" : "c",
          "coef" : 0.03,
          "link" : "cpd00001:-1"
        },          
        
        {
          "id"   : "cpd00129",
          "name" : "L-Proline",
          "comp" : "c",
          "coef" : 0.04,
          "link" : "cpd00001:-1"
        },  
        
        {
          "id"   : "cpd00054",
          "name" : "L-Serine",
          "comp" : "c",
          "coef" : 0.05,
          "link" : "cpd00001:-1"
        },          
        
        {
          "id"   : "cpd00161",
          "name" : "L-Threonine",
          "comp" : "c",
          "coef" : 0.06,
          "link" : "cpd00001:-1"
        },          
        
        {
          "id"   : "cpd00065",
          "name" : "L-Tryptophan",
          "comp" : "c",
          "coef" : 0.01,
          "link" : "cpd00001:-1"
        },          
        
        {
          "id"   : "cpd00069",
          "name" : "L-Tyrosine",
          "comp" : "c",
          "coef" : 0.02,
          "link" : "cpd00001:-1"
        },          
        
        {
          "id"   : "cpd00156",
          "name" : "L-Valine",
          "comp" : "c",
          "coef" : 0.08,
          "link" : "cpd00001:-1"
        }    
        
      ]
    },
    
    {
      "group_name" : "Inorganics",
      "mass" : 0.005,
      "unit_group" : "g",
      "unit_components" : "MOLSPLIT",
      "components" : [
        
        {
          "id"   : "cpd00048",
          "name" : "Sulfate",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd10515",
          "name" : "Fe2+",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd10516",
          "name" : "Fe3+",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00063",
          "name" : "Ca2+",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00058",
          "name" : "Cu2+",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00034",
          "name" : "Zn2+",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00254",
          "name" : "Mg+",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00205",
          "name" : "K+",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00099",
          "name" : "Cl-",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00030",
          "name" : "Mn2+",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00149",
          "name" : "Co2+",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00009",
          "name" : "Phosphate",
          "comp" : "c",
          "coef" : 1
        }
        
      ]
    },
    
    {
      "group_name" : "Cofactors",
      "mass" : 0.045,
      "unit_group" : "g",
      "unit_components" : "MOLSPLIT",
      "components" : [
        
        {
          "id"   : "cpd00042",
          "name" : "GSH",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00003",
          "name" : "NAD",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00017",
          "name" : "Adenosyl-L-methionine",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00010",
          "name" : "CoA",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00015",
          "name" : "FAD",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00056",
          "name" : "TPP",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd11493",
          "name" : "ACP",
          "comp" : "c",
          "coef" : 1,
          "link" : "cpd12370:-1"
        },
        
        {
          "id"   : "cpd00201",
          "name" : "10-Formyltetrahydrofolate",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00006",
          "name" : "NADP",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00345",
          "name" : "5-Methyltetrahydrofolate",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00016",
          "name" : "Pyridoxal-phosphate",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00087",
          "name" : "Tetrahydrofolate",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd00220",
          "name" : "Riboflavin",
          "comp" : "c",
          "coef" : 1
        }
        
      ]
    },
    
    {
      "group_name" : "Lipids",
      "mass" : 0.09,
      "unit_group" : "g",
      "unit_components" : "MOLSPLIT",
      "components" : [
        
        {
          "id"   : "cpd15540",
          "name" : "Phosphatidylglycerol-dioctadecanoyl",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd15793",
          "name" : "Stearoylcardiolipin-B-subtilis",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd15722",
          "name" : "Diisoheptadecanoylphosphatidylglycerol",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd15794",
          "name" : "Isoheptadecanoylcardiolipin-B-subtilis",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd15723",
          "name" : "Dianteisoheptadecanoylphosphatidylglycerol",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd15795",
          "name" : "Anteisoheptadecanoylcardiolipin-B-subtilis",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd15533",
          "name" : "phosphatidylethanolamine-dioctadecanoyl",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd15696",
          "name" : "Dianteisoheptadecanoylphosphatidylethanolamine",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd15695",
          "name" : "Diisoheptadecanoylphosphatidylethanolamine",
          "comp" : "c",
          "coef" : 1
        }
        
      ]
    },
    
    {
      "group_name" : "Cellwall",
      "mass" : 0.08,
      "unit_group" : "g",
      "unit_components" : "MOLSPLIT",
      "components" : [
        
        {
          "id"   : "cpd15432",
          "name" : "core-oligosaccharide-lipid-A",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd02229",
          "name" : "Bactoprenyl-diphosphate",
          "comp" : "c",
          "coef" : 1
        },
        
        {
          "id"   : "cpd15665",
          "name" : "Peptidoglycan-polymer-n-subunits",
          "comp" : "c",
          "coef" : 1,
          "link" : "cpd15666:-1"
        }
        
      ]
    }
  ]
}
```

*Explanation:* 

The basic fields for a biomass definition are:

- `id` - A a short ID for the biomass reaction
- `name` - Can be a longer name for the biomass reactions, that includes a brief description of the idea for the specific biomass.
- `ref` - Reference (if any).
- `engery_GAM` - How many ATP molecules are hydrolysed to ADP + Pi for growth-associated maintenance (GAM) per gram dry biomass formed? 
- `Domain` - Provide the domain for the organism type. (Bacteria/Archaea/Eukarya)
- `met_groups` - See below.

Biomass components should be assigned to groups (`met_groups`). Usually these groups are DNA, RNA, Protein, Others (e.g. inorganic components, lipids, cell wall components, etc...). The total mass of each group should be stated in the group parameter (`mass`) in the unit *gram*. The sum of all group masses should be 1 g.

The field `unit_components` can be "MOLFRACTION" or "MOLSPLIT". Both would result in the same effect: The molar coefficients for the group components are scaled to result in the exact mass in gram that is specified for the  corresponding group. The only difference between "MOLFRACTION" and "MOLSPLIT" is that  "MOLFRACTION" throws a warning if the summed coefficients of the group's components do not add up to the value of 1 (= 100%).

Then there are component entries. The format for each component is for instance:

```
        {
          "id"   : "cpd00161",
          "name" : "L-Threonine",
          "comp" : "c",
          "coef" : 0.0474316079511907,
          "link" : "cpd00001:-1"
        }, 
```

- `id` Refers to the metabolite ID in ModelSEED/gapseq
- `name` - Provide a name for the metabolite
- `c` - From which compartment should the metabolite be taken for the biomass? Usually this should be "c" for cytosol.
- `coef` - molar coefficient of the metabolite relative to all components within the same group
- `link` - When this metabolite is flowing into the biomass, should also other metabolites be consumed/produced with it?

The `link` is optional. For instance, if amino acids are elongated to peptides a water molecule is  produced with each amino acid added. The format in this case is `<compound_id>:<stoichiometry>`. In the above example: For every threonine incorporated into the biomass a water molecule (cpd00001) is released. In contrast, for nucleotides  that flow into DNA or RNA, pyrophospate (cpd00012) is produced with each nucleotide. Thus, the link here is `"link" : "cpd00012:-1"`.
