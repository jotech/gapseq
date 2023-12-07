# panDraft: species-level models

When dealing with environmental samples, it's common for Metagenome-Assembled Genomes (MAGs) obtained from genome centric metagenomics to be incomplete. Consequently, Genome-Scale Models (GEMs) derived from these incomplete MAGs lack a substantial portion of the metabolic potential of the corresponding species. In response, once a set of draft models is reconstructed, they can be conbined into a comprehensive model specific to the taxonomic group. This approach aims to fill the gaps present in individual models by combining homology-based searches and pan-reactome analysis.

Using `gapseq pan`, draft species-level models can be generated to enhance the representativeness of the species metabolism.

## Build the panDraft model of a species

Input need to be provided as **lists of pathways** to the desired files:
```
cat toy/panDraft/input4modelMerging.txt

toy/panDraft/MGYG000008797-draft.RDS
toy/panDraft/MGYG000009889-draft.RDS
toy/panDraft/MGYG000167774-draft.RDS
toy/panDraft/MGYG000173413-draft.RDS
```
The module requires the draft models of each MAG (.RDS), the reaction weight tables (-rxnWeigths.RDS), the gene to reaction associations tables (-rxnXgenes.RDS), and the pathways tables (-all-Pathways.tbl):
```
./gapseq pan -m toy/panDraft/input4modelMerging.txt -c toy/panDraft/input4WeightsMerging.txt -g toy/panDraft/input4XgenesMerging.txt -w toy/panDraft/input4Pathways.txt
```
The tool generates several outputs, including the draft model of the taxon, an updated version of the weight table, statistics on the species reactome, and `rxnXmod.tsv`, a binary matrix summarizing the presence/absence of reactions in each GEMs. Additionally, you can opt to perform a simple reactome comparison using the option `--only.binary.rxn.table TRUE`.

## Suggested number MAGs and mrf threshold value

We recommend utilizing a minimum of 30 MAGs for the reconstruction of a panDraft draft. However, there is no specified lower limit to this number.
Clearly, the minimum reaction frequency threshold (mrf), which determines whether a reaction should be included in the species-level model, is meaningful only when the number of MAGs used is above a certain minimum threshold. The parameter (mrf) can take values between 0 and 1 (default 0.07). A value of 0 means that all reactions present in any of the input models will be included in the draft model, while a value of 1 means that only reactions present in all input models will be included.

The mrf threshold can be modified using the option `--min.rxn.freq.in.mods`. 
```
./gapseq pan -m toy/panDraft/input4modelMerging.txt -c toy/panDraft/input4WeightsMerging.txt -g toy/panDraft/input4XgenesMerging.txt -w toy/panDraft/input4Pathways.txt --min.rxn.freq.in.mods 0.15
```