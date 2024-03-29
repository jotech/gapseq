# panDraft: species-level models

When dealing with environmental samples, it's common for Metagenome-Assembled Genomes (MAGs) obtained from genome centric metagenomics to be incomplete and contaminated. Consequently, Genome-Scale Models (GEMs) derived from these incomplete MAGs lack a substantial portion of the metabolic potential of the corresponding species. In response, once a set of draft models is reconstructed, they can be conbined into a comprehensive model specific to the taxonomic group. This approach aims to fill the gaps present in individual models by combining homology-based searches and pan-reactome analysis.

Using `gapseq pan`, draft species-level models can be generated to enhance the representativeness of the taxa (species) metabolism.

## Build the panDraft model of a species

The files required by *pan-Draft* include the draft models of the MAGs of interest (.RDS), the corresponding reaction weight tables (-rxnWeights.RDS), the gene-to-reaction association tables (-rxnXgenes.RDS), and the pathways tables (-all-Pathways.tbl). Here are two examples of how to run *pan-Draft* using four models for MAGs of *Faecalibacillus intestinalis*. These genomes were selected for demonstration purposes of the module's operation.

Input can be provided as:
1) lists of filepaths separated by commas
2) filepaths using wildcards
3) path to folders containing the desired files 
```
# Option 1:
./gapseq pan -m toy/MGYG000008797-draft.RDS,toy/MGYG000009889-draft.RDS,toy/MGYG000167774-draft.RDS,toy/MGYG000173413-draft.RDS -c toy/MGYG000008797-rxnWeights.RDS,toy/MGYG000009889-rxnWeights.RDS,toy/MGYG000167774-rxnWeights.RDS,toy/MGYG000173413-rxnWeights.RDS -g toy/MGYG000008797-rxnXgenes.RDS,toy/MGYG000009889-rxnXgenes.RDS,toy/MGYG000167774-rxnXgenes.RDS,toy/MGYG000173413-rxnXgenes.RDS -w toy/MGYG000008797-all-Pathways.tbl.gz,toy/MGYG000009889-all-Pathways.tbl.gz,toy/MGYG000167774-all-Pathways.tbl.gz,toy/MGYG000173413-all-Pathways.tbl.gz

# Option 2:
./gapseq pan -m toy/M*-draft.RDS -c toy/M*-rxnWeights.RDS -g toy/M*-rxnXgenes.RDS -w toy/M*-all-Pathways.tbl.gz

# Option 3:
mkdir toy/panDraft
mv toy/MGY* toy/panDraft
./gapseq pan -m toy/panDraft -c toy/panDraft -g toy/panDraft -w toy/panDraft
```
The tool generates several outputs, including the draft model of the taxon (`panModel-draft.RDS`), an updated version of the weights (`panModel-rxnWeigths.RDS`), gene-to-reaction (`panModel-rxnXgenes.RDS`) and pathway table (`panModel-tmp-Pathways.tbl`), statistics on the species reactome features (`pan-reactome_stat.tsv`), and a binary matrix summarizing the presence/absence of reactions in each GEMs (`rxnXmod.tsv`). 

Additionally, you can opt to perform a simple reactome comparison using the option `--only.binary.rxn.tbl`. This will generate only the binary matrix summarizing the presence/absence of reactions.
```
./gapseq pan -m toy/M*-draft.RDS --only.binary.rxn.tbl TRUE
```
### Suggested number MAGs and mrf threshold value

We recommend utilizing a minimum of 30 MAGs for the reconstruction of a panDraft draft. However, there is no specified lower limit to this number.
Clearly, the minimum reaction frequency threshold (mrf), which determines whether a reaction should be included in the species-level model, is meaningful only when the number of MAGs used is above a certain minimum number. The parameter (mrf) can take values between 0 and 1 (default 0.07). A value of 0 means that all reactions present in any of the input models will be included in the draft model, while a value of 1 means that only reactions present in all input models will be included.

The mrf threshold can be modified using the option `--min.rxn.freq.in.mods`. 
```
./gapseq pan -m toy/M*-draft.RDS -c toy/M*-rxnWeights.RDS -g toy/M*-rxnXgenes.RDS -w toy/M*-all-Pathways.tbl.gz --min.rxn.freq.in.mods 0.15
```

## Integration with further features of gapseq

### Generate a functional model

After reconstructing the draft species-level models of the taxa of interest, the gapseq pipeline can proceed with the gap filling step. In this case, the input to the `gapseq fill` module should be the updated version of the reaction weights and reactions to gene tables.
```
./gapseq fill -m toy/panModel-draft.RDS -c toy/panModel-rxnWeigths.RDS -g toy/panModel-rxnXgenes.RDS -n dat/media/TSBmed.csv
```
The resulting model is functional and can be used for further metabolic model simulations.

### Predict the growing medium

The output of **pan-Draft** can be also integrated with the module `gapseq medium` to predict the growing medium of the species-level model. Here it is an example:
```
# Using the species-level draft
./gapseq medium -m panModel-draft.RDS -p panModel-tmp-Pathways.tbl

# Using the species-level functional model
./gapseq medium -m panModel.RDS -p panModel-tmp-Pathways.tbl
```
