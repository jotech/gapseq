# Reconstruction of metabolic networks for archaea 

> <mark>Please note</mark>: This tutorial assumes a gapseq  <= v1.3.1. An updated version of the tutorial for gapseq >= v1.4.0 is under construction.

Archaea are prokaryotes with important unique characteristics which separate them from the other domains of life.
Certain environmental process such as the production of methane is only described for archaea and there are estimates that they could make up to more than 10% of the human gut microbiome [1](https://doi.org/10.1016/j.anaerobe.2011.03.001).
This tutorial is about the reconstruction and analysis of methane-producing archaea.

## Methanosarcina barkeri
`Methanosarcina barkeri` is a methane-producing archaea that lives in sewage, mud, or rumen and is able to catabolize a variety of carbon sources.
A genome of `M. barkeri` is available on ncbi:
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/895/GCF_000195895.1_ASM19589v1/GCF_000195895.1_ASM19589v1_genomic.fna.gz
```

The following gapseq commands will search for pathways, transporters and will create a metabolic model for `M. barkeri`:
```
./gapseq find -p all -t Archaea GCF_000195895.1_ASM19589v1_genomic.fna.gz
./gapseq find-transport GCF_000195895.1_ASM19589v1_genomic.fna.gz
./gapseq draft -r GCF_000195895.1_ASM19589v1_genomic-all-Reactions.tbl -t GCF_000195895.1_ASM19589v1_genomic-Transporter.tbl -b archaea -u 200 -l 100 -a 1 -p GCF_000195895.1_ASM19589v1_genomic-all-Pathways.tbl
./gapseq fill -m GCF_000195895.1_ASM19589v1_genomic-draft.RDS -n ./dat/media/MM_anaerobic_CO2_H2.csv -c GCF_000195895.1_ASM19589v1_genomic-rxnWeights.RDS -b 100 -g GCF_000195895.1_ASM19589v1_genomic-rxnXgenes.RDS
```
The metabolic model is created as SBML or RDS file and can further be used for constraint-based analysis. The following example is using R.
```R
library(sybil)
library(data.table)
mod <- readRDS("GCF_000195895.1_ASM19589v1_genomic.RDS")
mtf <- pfba(mod); print(mtf)
ex <- findExchReact(mod)
ex.dt <- data.table(ex=ex@react_id, met=mod@met_name[ex@met_pos], flux=mtf@fluxdist@fluxes[ex@react_pos])
ex.dt[flux<0][order(flux)] # substrate: h2+co2
ex.dt[flux>0][order(flux)] # products: h2o+ch4
```
As you can see, the metabolic model of `M. barkeri` predicts, beside an growth rate of 0.15/h, the uptake of H2 and CO2 to produce methane (CH4) and water . 

## Methanogens with and without cytochromes
An interesting criterion for distinguishing methanogens is the presence of cytochromes [2](http://dx.doi.org/10.1038/nrmicro1931).
Methanogens with cytochromes use methanophenanzine (a menaquinone analogue) to build up a membrane potential for energy conservation via electron transport phosphorylation.
Methanogens without cytochromes on the other hand do not use methanophenanzine.
Instead, they rely on flavin-based electron bifurcation for energy conservation.
`M. barkeri` is a methanogen with cytochromes and the reactions that used 'methanophenanzine carry active fluxes:
```R
sol.dt <- data.table(rxn=mod@react_id, name=mod@react_name, flux=mtf@fluxdist@fluxes[1:mod@react_num])
sol.dt[grepl("rxn90059|rxn90058", rxn)] # methanophenanzine-using reactions
```

Let's consider another methanogenic archaeon which does not possess cytochromes.
Methanothermobacter thermautotrophicus is a thermophilic methanogen found in sewage sludge [3](https://microbewiki.kenyon.edu/index.php/Methanothermobacter_thermautotrophicus).
The model creation with gapseq is just like before: 
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/645/GCF_000008645.1_ASM864v1/GCF_000008645.1_ASM864v1_genomic.fna.gz
./gapseq find -p all -t Archaea GCF_000008645.1_ASM864v1_genomic.fna.gz
./gapseq find-transport GCF_000008645.1_ASM864v1_genomic.fna.gz
./gapseq draft -r GCF_000008645.1_ASM864v1_genomic-all-Reactions.tbl -t GCF_000008645.1_ASM864v1_genomic-Transporter.tbl -b archaea -u 200 -l 100 -a 1 -p GCF_000008645.1_ASM864v1_genomic-all-Pathways.tbl
./gapseq fill -m GCF_000008645.1_ASM864v1_genomic-draft.RDS -n dat/media/MM_anaerobic_CO2_H2.csv -c GCF_000008645.1_ASM864v1_genomic-rxnWeights.RDS -b 100 -g GCF_000008645.1_ASM864v1_genomic-rxnXgenes.RDS
```
The prediction of substrates and products for `M. thermautotrophicus` is the same as for `M. barkeri`:
```R
mod2 <- readRDS("GCF_000008645.1_ASM864v1_genomic.RDS")
mtf2 <- optimizeProb(mod2, algorithm="mtf"); print(mtf2)
ex2 <- findExchReact(mod2)
ex2.dt <- data.table(ex=ex2@react_id, met=mod2@met_name[ex2@met_pos], flux=mtf2@fluxdist@fluxes[ex2@react_pos])
ex2.dt[flux<0][order(flux)] # substrate: h2+co2
ex2.dt[flux>0][order(flux)] # products: h2o+ch4
```
Nonetheless, the underlying energy conserving process is completely different.
The model of `M. thermautotrophicus` uses an electron bifurcating reaction that couples the energetic unfavorable electron transfer from H2 to ferredoxin with the energetic favorable reduction of CoM-S-S-CoB:
```R
sol2.dt <- data.table(rxn=mod2@react_id, name=mod2@react_name, flux=mtf2@fluxdist@fluxes[1:mod2@react_num])
sol2.dt[grepl("rxn42401", rxn)] # electron bifurcating reaction
```

Interestingly, the methane yield of methanogens with cytochromes is higher than the methane yield of methanogens without cytochromes.
This can be reproduced by the metabolic models:
```R
mtf@fluxdist@fluxes[which(mod@obj_coef!=0)] / mtf@fluxdist@fluxes[grep("EX_cpd01024", mod@react_id)] * 1000 # yield_ch4 = 6.91 / mol
mtf2@fluxdist@fluxes[which(mod2@obj_coef!=0)] / mtf2@fluxdist@fluxes[grep("EX_cpd01024", mod2@react_id)] * 1000 # yield_ch4 = 3.33 / mol
```
The methane yield of both organisms are close to reported experimental values [4](http://dx.doi.org/10.1038/nrmicro1931).
