# Lactic acid bacteria in yogurt production

*Lactobacillus delbrueckii* subsp. *bulgaricus* and *Streptococcus thermophilus* are both lactic acid bacteria, which are frequently used for the production of yogurt. The organisms differ in their lactose degradation and their fermentation end products. In this tutorial, genome-scale metabolic models will be reconstructed using **gapseq** starting from the organisms' genomes in multi-protein sequences fasta file (translated sequences of protein-coding genes).

##### Input files

- Genome assemblies:
  - *Lactobacillus delbrueckii* subsp. *bulgaricus* ATCC 11842 = JCM 1002. RefSeq: `GCF_000056065.1`
  - *Streptococcus thermophilus* ATCC 19258. RefSeq: `GCF_010120595.1`
- Growth media (milk) file: `milk.csv` 
  This growth media is based on the main ingredients described for whole milk (without fatty acids), as listed [here](https://frida.fooddata.dk/food/1265?lang=en).

##### Preparations

Download required files and rename assembly files for the ease of this tutorial.

```sh
mkdir yoghurt
cd yoghurt

# Download genome assemblies from NCBI (RefSeq) in protein sequences
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/056/065/GCF_000056065.1_ASM5606v1/GCF_000056065.1_ASM5606v1_protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/010/120/595/GCF_010120595.1_ASM1012059v1/GCF_010120595.1_ASM1012059v1_protein.faa.gz

# Download gapfill-medium file from github
wget https://github.com/Waschina/gapseq.tutorial.data/raw/master/yogurt/milk.csv

# Rename genomes
mv GCF_000056065.1_ASM5606v1_protein.faa.gz ldel.faa.gz
mv GCF_010120595.1_ASM1012059v1_protein.faa.gz sthe.faa.gz
```

##### gapseq reconstruction pipeline

(1) Reaction & pathway prediction
(2) Transporter prediction
(3) Draft model reconstruction
(4) Gapfilling

```sh
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

```R
getMetaboliteProduction <- function(mod) {
  require(cobrar)
  require(data.table)
  
  # MTF
  sol.mtf <- pfba(mod)
  dt.mtf  <- data.table(ex = mod@react_id,
                        name = mod@react_name,
                        mtf.flux = sol.mtf@fluxes)
  
  # FVA
  sol.fv <- fva(mod, react = mod@react_id[grep("^EX_cpd[0-9]+_e0", mod@react_id)])
  sol.fv$growth.fraction <- NULL
  
  dt <- merge(dt.mtf, sol.fv, by.x = "ex", by.y = "react")
  dt <- dt[mtf.flux > 1e-6]
  dt <- dt[, .(ex, rxn.name = name, l = min.flux, u = max.flux, mtf.flux)]
  
  return(dt[order(-mtf.flux)])
}
```

Now, we can apply this function to the network models of *L. delbrueckii* and *S. thermophilus* to predict the top 10 produced metabolic by-products.

```R
library(cobrar)

ld <- readRDS("ldel.RDS") # for L. delbrueckii
st <- readRDS("sthe.RDS") # for S. thermophilus

getMetaboliteProduction(ld)[1:10]
getMetaboliteProduction(st)[1:10]
```

Output for *L. delbrueckii* (`ld`):

```
                ex               rxn.name          l          u   mtf.flux
            <char>                 <char>      <num>      <num>      <num>
 1: EX_cpd00067_e0         H+-e0 Exchange 4.63973657 4.74111902 4.63973657
 2: EX_cpd00159_e0  L-Lactate-e0 Exchange 0.00000000 4.60791680 4.53569122
 3: EX_cpd00108_e0  Galactose-e0 Exchange 2.50000000 2.50000000 2.50000000
 4: EX_cpd00011_e0        CO2-e0 Exchange 0.25028483 0.32251041 0.32251041
 5: EX_cpd00141_e0 Propionate-e0 Exchange 0.16796159 0.21865281 0.16796159
 6: EX_cpd00047_e0    Formate-e0 Exchange 0.16354480 0.21423602 0.16354480
 7: EX_cpd00239_e0        H2S-e0 Exchange 0.09359832 0.09359832 0.09359832
 8: EX_cpd00324_e0       MTTL-e0 Exchange 0.08643477 0.08643477 0.08643477
 9: EX_cpd00013_e0        NH3-e0 Exchange 0.08472240 0.18610486 0.08472240
10: EX_cpd00029_e0    Acetate-e0 Exchange 0.00000000 0.02153435 0.02153435

```

Output for S. thermophilus (`st`):

```
                ex                    rxn.name            l           u   mtf.flux
            <char>                      <char>        <num>       <num>      <num>
 1: EX_cpd00067_e0              H+-e0 Exchange  0.265336716 18.24558254 9.54662000
 2: EX_cpd00159_e0       L-Lactate-e0 Exchange  0.000000000  9.03085323 8.30919630
 3: EX_cpd00047_e0         Formate-e0 Exchange  0.184869492  7.53849623 0.88018350
 4: EX_cpd00029_e0         Acetate-e0 Exchange  0.000000000  9.37480074 0.43403457
 5: EX_cpd00141_e0      Propionate-e0 Exchange  0.077254201  2.58930860 0.12591620
 6: EX_cpd00011_e0             CO2-e0 Exchange  0.000000000  9.39714849 0.08666366
 7: EX_cpd00239_e0             H2S-e0 Exchange  0.000000000  0.08858863 0.08643612
 8: EX_cpd00324_e0            MTTL-e0 Exchange  0.075101685  0.16369032 0.07725420
 9: EX_cpd00066_e0 L-Phenylalanine-e0 Exchange -0.033081729  0.02681010 0.02681010
10: EX_cpd00036_e0       Succinate-e0 Exchange  0.003587527  0.47281948 0.02392311
```

As expected, both organisms produce Lactate in the FBA+MTF solution. In contrast to *S. thermophilus*, the FBA simulation predicted a release of Galactose by *L. debrueckii*. In fact, *L. debrueckii* is usually reported to be Galactose negative; i.e. does not produce acid from this hexose ([https://bacdive.dsmz.de/strain/6449](https://bacdive.dsmz.de/strain/6449)) and utilized only the glucose part of lactose, while *S. thermophilus* has been reported to be Galactose positive ([https://bacdive.dsmz.de/strain/14786](https://bacdive.dsmz.de/strain/14786)).
