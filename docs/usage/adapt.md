# Adapting models

Once a model is reconstructed and gap filled it is ready to predict cellular growth for example using flux balance analysis.
But still there can be shortcomings and so that experimental phenotypes are not met by simulations.
By `gapseq adapt`, models can be manually curated and extended to improve the overall accuracy of an metabolic model.

## Change reactions or entire pathways of a model
Reactions can be added by multiple identifiers:
```
./gapseq adapt -m toy/myb71.RDS -a rxn00814

# KEGG
./gapseq adapt -m toy/myb71.RDS -a R01098

# MetaCYC
./gapseq adapt -m toy/myb71.RDS -a NAPHTHALENE-12-DIOXYGENASE-RXN
Loading model files toy/myb71.RDS
```

In the same way, reactions can be removed:
```
./gapseq adapt -m toy/myb71.RDS -r rxn00001

```

In addition, reactions of entire pathways can be integrated.
For example, the photosynthesis pathway could be integrated into E. coli:
```
./gapseq adapt -a PHOTOALL-PWY -m toy/ecoli.RDS
Loading model files toy/ecoli.RDS
```
Or removed
```
./gapseq adapt -m toy/ecoli.RDS -r GLYCOLYSIS
```

## Modify growth conditions of a model
Phenotypic growth data of an organism can be available for example from literature or from screening experiments like Biolog Microplates.
If reconstructed and gap filled models do not match the known growth behavior, gapseq is able to try fixing the model.
This is be done by either adding reactions to allow growth under one known specific growth condition or by removing reactions that permit model growth under conditions in which there should be no growth possible.

By using `gapseq adapt`, a substance can be defined which a model should be able to use. To enable butyrate (cpd00211) uptake for example try the following:
```
./gapseq adapt -m toy/myb71.RDS -w cpd00211:TRUE -c toy/myb71-rxnWeights.RDS -g toy/myb71-rxnXgenes.RDS -b toy/myb71-all-Reactions.tbl
```

Similarly, the ability to metabolize a substance can also be removed from a model. To remove the usage of glucose use:
```
./gapseq adapt -m toy/myb71.RDS -w cpd00027:FALSE -c toy/myb71-rxnWeights.RDS -g toy/myb71-rxnXgenes.RDS -b toy/myb71-all-Reactions.tbl
```

Both, enabling and disabling of growth phenotypes can be combined in a single command providing a list of substances with there growth state:
```
./gapseq adapt -m toy/myb71.RDS -w cpd00211:TRUE,cpd00027:FALSE -c toy/myb71-rxnWeights.RDS -g toy/myb71-rxnXgenes.RDS -b toy/myb71-all-Reactions.tbl
```
