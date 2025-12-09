# Lactic acid bacteria in yogurt production
2025-12-04

*Lactobacillus delbrueckii* subsp. *bulgaricus* and *Streptococcus thermophilus* are both lactic acid bacteria, which are frequently used for the production of yogurt. The organisms differ in their lactose degradation and their fermentation end products. In this tutorial, genome-scale metabolic models will be reconstructed using **gapseq** starting from the organisms' genomes in multi-protein sequences fasta file (translated sequences of protein-coding genes).

##### Input files

- Genome assemblies:
  - *Lactobacillus delbrueckii* subsp. *bulgaricus* ATCC 11842 = JCM 1002. RefSeq: `GCF_000056065.1`
  - *Streptococcus thermophilus* ATCC 19258. RefSeq: `GCF_010120595.1`
- Growth media (milk) file: `milk.csv` 
  This growth media is based on the main ingredients described for whole milk (without fatty acids), as listed [here](https://frida.fooddata.dk/food/1265?lang=en).

##### Preparations

Download required files and rename assembly files for the ease of this tutorial.


``` sh
# Download genome assemblies from NCBI (RefSeq) in protein sequences
wget -q -O ldel.faa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/056/065/GCF_000056065.1_ASM5606v1/GCF_000056065.1_ASM5606v1_protein.faa.gz
wget -q -O sthe.faa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/010/120/595/GCF_010120595.1_ASM1012059v1/GCF_010120595.1_ASM1012059v1_protein.faa.gz

# Download gapfill-medium file from github
wget -q https://github.com/Waschina/gapseq.tutorial.data/raw/master/yogurt/milk.csv
```

##### gapseq reconstruction pipeline

(1) Reaction & pathway prediction
(2) Transporter prediction
(3) Draft model reconstruction
(4) Gapfilling


``` sh
modelA="ldel"
modelB="sthe"

# (1) Reaction & Pathway prediction
gapseq find -p all -b 200 -m auto -t auto -A diamond $modelA.faa.gz
gapseq find -p all -b 200 -m auto -t auto -A diamond $modelB.faa.gz

# (2) Transporter prediction
gapseq find-transport -b 200 -A diamond $modelA.faa.gz 
gapseq find-transport -b 200 -A diamond $modelB.faa.gz

# (3) Building Draft Model - based on Reaction-, Pathway-, and Transporter prediction
gapseq draft -r $modelA-all-Reactions.tbl -t $modelA-Transporter.tbl -p $modelA-all-Pathways.tbl -u 200 -l 100
gapseq draft -r $modelB-all-Reactions.tbl -t $modelB-Transporter.tbl -p $modelB-all-Pathways.tbl -u 200 -l 100

# (4) Gapfilling
gapseq fill -m $modelA-draft.RDS -n milk.csv -b 100
gapseq fill -m $modelB-draft.RDS -n milk.csv -b 100
```

##### FBA and FVA prediction of metabolic by-products

Here, we will use the R-Package [`cobrar`](https://waschina.github.io/cobrar/) to perform Flux-Balance-Analysis (FBA) and Flux-Variability-Analysis (FVA) with the two reconstructed network models.

First, we define a function, that automatically performs FBA with the minimization of total flux (MTF) as secondary objective and FVA for all exchange reactions. The function also summarizes the results in a sorted `data.table`.


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

getMetaboliteProduction <- function(mod) {
  sol.mtf <- pfba(mod)
  dt.mtf  <- data.table(getExchanges(mod, sol.mtf))
  dt.mtf <- dt.mtf[flux >= 1e-5 & ID != "EX_cpd11416_c0"]
  
  return(dt.mtf[order(-flux)])
}
```

Now, we can apply this function to the network models of *L. delbrueckii* and *S. thermophilus* to predict the top 10 produced metabolic by-products.

Results for *L. delbrueckii* (`ld`):


``` r
ld <- readRDS("ldel.RDS") # for L. delbrueckii
getMetaboliteProduction(ld)[1:10]
```

```
##                 ID                   name       flux
##             <char>                 <char>      <num>
##  1: EX_cpd00067_e0         H+-e0 Exchange 4.65237700
##  2: EX_cpd00159_e0  L-Lactate-e0 Exchange 4.55008609
##  3: EX_cpd00108_e0  Galactose-e0 Exchange 2.50000000
##  4: EX_cpd00011_e0        CO2-e0 Exchange 0.30908772
##  5: EX_cpd00141_e0 Propionate-e0 Exchange 0.16793024
##  6: EX_cpd00047_e0    Formate-e0 Exchange 0.16350913
##  7: EX_cpd00013_e0        NH3-e0 Exchange 0.10052070
##  8: EX_cpd00239_e0        H2S-e0 Exchange 0.09359205
##  9: EX_cpd00324_e0       MTTL-e0 Exchange 0.08642149
## 10: EX_cpd00130_e0   L-Malate-e0 Exchange 0.02223483
```


Results for S. thermophilus (`st`):


``` r
st <- readRDS("sthe.RDS") # for S. thermophilus
getMetaboliteProduction(st)[1:10]
```

```
##                 ID                        name       flux
##             <char>                      <char>      <num>
##  1: EX_cpd00067_e0              H+-e0 Exchange 9.53844835
##  2: EX_cpd00159_e0       L-Lactate-e0 Exchange 8.30267113
##  3: EX_cpd00047_e0         Formate-e0 Exchange 0.87912596
##  4: EX_cpd00029_e0         Acetate-e0 Exchange 0.43391447
##  5: EX_cpd00141_e0      Propionate-e0 Exchange 0.12535533
##  6: EX_cpd00239_e0             H2S-e0 Exchange 0.08644827
##  7: EX_cpd00011_e0             CO2-e0 Exchange 0.08595875
##  8: EX_cpd00324_e0            MTTL-e0 Exchange 0.07727458
##  9: EX_cpd00066_e0 L-Phenylalanine-e0 Exchange 0.02678608
## 10: EX_cpd00036_e0       Succinate-e0 Exchange 0.02390168
```


As expected, both organisms produce lactate in the FBA+MTF solution. In contrast to *S. thermophilus*, the FBA simulation predicted a release of galactose by *L. debrueckii*. In fact, *L. debrueckii* is usually reported to be galactose-negative; i.e. does not produce acid from this hexose ([https://bacdive.dsmz.de/strain/6449](https://bacdive.dsmz.de/strain/6449)) and utilized only the glucose part of lactose, while *S. thermophilus* has been reported to be galactose-positive ([https://bacdive.dsmz.de/strain/14786](https://bacdive.dsmz.de/strain/14786)).


