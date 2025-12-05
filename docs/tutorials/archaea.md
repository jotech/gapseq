# Reconstruction of metabolic networks for archaea 
2025-12-04

Archaea are prokaryotes with important unique characteristics which separate them from the other domains of life.
Certain environmental process such as the production of methane is only described for archaea and there are estimates that they could make up to more than 10% of the human gut microbiome[[1](https://doi.org/10.1016/j.anaerobe.2011.03.001)].

This tutorial is about the reconstruction and analysis of methane-producing archaea.

##### *Methanosarcina barkeri*

*Methanosarcina barkeri* is a methane-producing archaea that lives in sewage, mud, and rumen and is able to catabolize a variety of carbon sources.
A genome of *M. barkeri* is available on NCBI:


``` sh
wget -q -O mbar.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/970/025/GCF_000970025.1_ASM97002v1/GCF_000970025.1_ASM97002v1_genomic.fna.gz
```

The following gapseq commands will search for pathways, transporters and will create a metabolic model for *M. barkeri*:


``` sh
gapseq find -p all -t Archaea -m Archaea -A diamond mbar.fna.gz
gapseq find-transport -A diamond mbar.fna.gz
gapseq draft -r mbar-all-Reactions.tbl -t mbar-Transporter.tbl -u 200 -l 100 -p mbar-all-Pathways.tbl
gapseq fill -m mbar-draft.RDS -n ../../dat/media/MM_anaerobic_CO2_H2.csv -b 100 -e highH2
```

The metabolic model is created as SBML and RDS file and can further be used for constraint-based analysis. The following example is using R and the extension package cobrar.


``` r
library(cobrar)
```

```
## Loading required package: Matrix
```

```
## cobrar uses...
##  - libSBML (v. 5.20.4)
##  - glpk (v. 5.0)
```

``` r
library(data.table)

mod <- readRDS("mbar.RDS")
mtf <- pfba(mod); print(mtf)
```

```
## Algorithm:              pFBA 
## Solver status:          solution is optimal 
## Optimization status:    optimization process was successful 
## Objective fct. value:   0.3560925 
## Secondary objective:    811.9304
```

``` r
ex.dt <- data.table(getExchanges(mod, mtf))
print("Substrates:"); ex.dt[flux<0][order(flux)] # substrate: h2+co2
```

```
## [1] "Substrates:"
```

```
##                ID                  name          flux
##            <char>                <char>         <num>
## 1: EX_cpd11640_e0        H2-e0 Exchange -1.000000e+02
## 2: EX_cpd00011_e0       CO2-e0 Exchange -3.327725e+01
## 3: EX_cpd00013_e0       NH3-e0 Exchange -4.150228e+00
## 4: EX_cpd00009_e0 Phosphate-e0 Exchange -3.464049e-01
## 5: EX_cpd00048_e0   Sulfate-e0 Exchange -9.797477e-02
## 6: EX_cpd00149_e0      Co2+-e0 Exchange -1.570004e-03
## 7: EX_cpd00244_e0      Ni2+-e0 Exchange -6.680868e-04
```

``` r
print("Products:"); ex.dt[flux>0][order(-flux)] # products: h2o+ch4
```

```
## [1] "Products:"
```

```
##                ID                     name         flux
##            <char>                   <char>        <num>
## 1: EX_cpd00001_e0          H2O-e0 Exchange 6.068646e+01
## 2: EX_cpd01024_e0      Methane-e0 Exchange 1.659132e+01
## 3: EX_cpd00067_e0           H+-e0 Exchange 4.124422e+00
## 4: EX_cpd00036_e0    Succinate-e0 Exchange 4.349951e-01
## 5: EX_cpd00055_e0 Formaldehyde-e0 Exchange 4.274185e-01
## 6: EX_cpd11416_c0           EX cpd11416 c0 3.560925e-01
## 7: EX_cpd00073_e0         Urea-e0 Exchange 1.189194e-02
## 8: EX_cpd00180_e0      Oxalate-e0 Exchange 2.839369e-04
```
As you can see, the metabolic model of *M. barkeri* predicts, beside a growth rate of 0.36/h, the uptake of H<sub>2</sub> (-100 mmol/gDW/h) and CO<sub>2</sub> (-33 mmol/gDW/h) to produce methane (17 mmol/gDW/h) and water (61 mmol/gDW/h).

##### Methanogens with and without cytochromes

An interesting criterion for distinguishing methanogens is the presence of cytochromes[[2](http://dx.doi.org/10.1038/nrmicro1931)].
Methanogens with cytochromes use methanophenanzine (a menaquinone analogue) to build up a membrane potential for energy conservation via electron transport phosphorylation.

Methanogens without cytochromes on the other hand do not use methanophenanzine. Instead, they rely on flavin-based electron bifurcation for energy conservation.

*M. barkeri* is a methanogen with cytochromes and the reactions that used methanophenanzine carry active fluxes:




``` r
sol.dt <- data.table(rxn=mod@react_id, name=mod@react_name, flux=mtf@fluxes, equation = printReaction(mod, mod@react_id))
sol.dt[grepl("rxn90059|rxn90058", rxn)] # methanophenanzine-using reactions
```

```
##            rxn                                                                           name     flux                                                                                                                           equation
##         <char>                                                                         <char>    <num>                                                                                                                             <char>
## 1: rxn90058_c0           Coenzyme B:coenzyme M:methanophenazine oxidoreductase (chemiosmosis) 16.59132 (3) H+-c0 + (1) CoM-S-S-CoB-c0 + (1) Dihydromethanophenazine-c0 <==> (1) CoM-c0 + (1) HTP-c0 + (1) Methanophenazine-c0 + (2) H+-e0
## 2: rxn90059_c0 hydrogen:2-(2,3-dihydropentaprenyloxy)phenazine oxidoreductase  (chemiosmosis) 16.59132                                     (2) H+-c0 + (1) Methanophenazine-c0 + (1) H2-c0 --> (1) Dihydromethanophenazine-c0 + (2) H+-e0
```



Let's consider another methanogenic archaeon which does not possess cytochromes. *Methanothermobacter thermautotrophicus* is a thermophilic methanogen found in sewage sludge[[3](https://microbewiki.kenyon.edu/index.php/Methanothermobacter_thermautotrophicus)].

The model reconstruction with gapseq is just like before: 


``` sh
wget -q -O mthe.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/645/GCF_000008645.1_ASM864v1/GCF_000008645.1_ASM864v1_genomic.fna.gz

gapseq find -p all -t Archaea -m Archaea -A diamond mthe.fna.gz
gapseq find-transport -A diamond mthe.fna.gz
gapseq draft -r mthe-all-Reactions.tbl -t mthe-Transporter.tbl -u 200 -l 100 -p mthe-all-Pathways.tbl
gapseq fill -m mthe-draft.RDS -n ../../dat/media/MM_anaerobic_CO2_H2.csv -b 100 -e highH2
```

The prediction of substrates and products for *M. thermautotrophicus* is the same as for *M. barkeri*:


``` r
mod2 <- readRDS("mthe.RDS")
mtf2 <- pfba(mod2); print(mtf2)
```

```
## Algorithm:              pFBA 
## Solver status:          solution is optimal 
## Optimization status:    optimization process was successful 
## Objective fct. value:   0.1225659 
## Secondary objective:    747.8264
```

``` r
ex2.dt <- data.table(getExchanges(mod2, mtf2))
print("Substrates:"); ex2.dt[flux<0][order(flux)] # substrate: h2+co2
```

```
## [1] "Substrates:"
```

```
##                ID                  name          flux
##            <char>                <char>         <num>
## 1: EX_cpd11640_e0        H2-e0 Exchange -1.000000e+02
## 2: EX_cpd00011_e0       CO2-e0 Exchange -2.755384e+01
## 3: EX_cpd00013_e0       NH3-e0 Exchange -1.428495e+00
## 4: EX_cpd00009_e0 Phosphate-e0 Exchange -1.192315e-01
## 5: EX_cpd00048_e0      Sulfate Exchange -3.372261e-02
## 6: EX_cpd00149_e0      Co2+-e0 Exchange -5.403904e-04
## 7: EX_cpd00244_e0      Ni2+-e0 Exchange -2.299534e-04
```

``` r
print("Products:"); ex2.dt[flux>0][order(-flux)] # products: h2o+ch4
```

```
## [1] "Products:"
```

```
##                ID                     name         flux
##            <char>                   <char>        <num>
## 1: EX_cpd00001_e0          H2O-e0 Exchange 5.361540e+01
## 2: EX_cpd01024_e0      Methane-e0 Exchange 2.233807e+01
## 3: EX_cpd00067_e0           H+-e0 Exchange 1.104639e+00
## 4: EX_cpd00055_e0 Formaldehyde-e0 Exchange 1.775020e-01
## 5: EX_cpd11416_c0           EX cpd11416 c0 1.225659e-01
## 6: EX_cpd00204_e0           CO-e0 Exchange 2.402992e-02
## 7: EX_cpd00029_e0      Acetate-e0 Exchange 8.503555e-03
## 8: EX_cpd00073_e0         Urea-e0 Exchange 4.093170e-03
## 9: EX_cpd00180_e0      Oxalate-e0 Exchange 9.773018e-05
```

Altough the main substrates and metabolic products are the same, the underlying energy conserving process is completely different. The model of *M. thermautotrophicus* uses an electron bifurcating reaction that couples the energetic unfavorable electron transfer from H<sub>2</sub> to ferredoxin with the energetic favorable reduction of CoM-S-S-CoB:


``` r
sol2.dt <- data.table(rxn=mod2@react_id, name=mod2@react_name, flux=mtf2@fluxes, equation = printReaction(mod2, mod2@react_id))
sol2.dt[grepl("rxn42401", rxn)] # electron bifurcating reaction
```

```
##            rxn                                            name      flux
##         <char>                                          <char>     <num>
## 1: rxn42401_c0 H2:CoB-CoM heterodisulfide,ferredoxin reductase -22.33807
##                                                                                                                         equation
##                                                                                                                           <char>
## 1: (1) H+-c0 + (1) CoM-c0 + (1) HTP-c0 + (2) Reducedferredoxin-c0 <-- (1) CoM-S-S-CoB-c0 + (2) Oxidizedferredoxin-c0 + (2) H2-c0
```

Interestingly, the growth yield $Y_{CH_4}$ (i.e, gramm dry weight (gDW) per mole of produced methane) of methanogens with cytochromes are is reported to be higher than the growth yield of methanogens without cytochromes[[4](http://dx.doi.org/10.1038/nrmicro1931)].

This can be reproduced by the metabolic models. Let's calculate the growth yields $Y_{CH_4}$ for both methanogens:


``` r
yield_mbar <- mtf@obj / mtf@fluxes[react_pos(mod,"EX_cpd01024_e0")] * 1000
yield_mbar
```

```
## [1] 21.46258
```

``` r
yield_mthe <- mtf2@obj / mtf2@fluxes[react_pos(mod2,"EX_cpd01024_e0")] * 1000
yield_mthe
```

```
## [1] 5.486862
```

The models predict, that the growth yield of *M. barkeri* is with 21.5 gDW/mol about 3.9 times higher than the growth yield of *M. thermautotrophicus* (5.5 gDW/mol), which is close to the ratio observed in experimental studies[[4](http://dx.doi.org/10.1038/nrmicro1931)].





