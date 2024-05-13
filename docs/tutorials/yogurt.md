# Lactic acid bacteria in yogurt production

*Lactobacillus delbrueckii* subsp. *bulgaricus* and *Streptococcus thermophilus* are both lactic acid bacteria, which are frequently used for the production of yogurt. The organisms differ in their lactose degradation and their fermentation end products. In this tutorial, genome-scale metabolic models will be reconstructed using **gapseq** starting from the organisms' genomes in multi-protein sequences fasta file (translated sequences of protein-coding genes).

Required gapseq version: 1.2 0d0123e (or later) 

##### Input files

- Genome assemblies:
  - *Lactobacillus delbrueckii* subsp. *bulgaricus* ATCC 11842 = JCM 1002. RefSeq: `GCF_000056065.1`
  - *Streptococcus thermophilus* ATCC 19258. RefSeq: `GCF_010120595.1`
- Growth media (milk) file: `milk.csv` 
  This growth media is based on the main ingredients described for whole milk (without fatty acids), as listed [here](https://frida.fooddata.dk/food/1265?lang=en).

>  NOTE: All intermediate files produced by the commands below are stored at the github repository [gapseq.tutorial.data](https://github.com/Waschina/gapseq.tutorial.data), which you can download/clone if you wish to start not at the beginning but at a later step of this tutorial or want to cross-check your results.



##### Preparations

Download required files and rename assembly filed for the ease of this tutorial.

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
gapseq find -p all -b 200 -m auto -t auto $modelA.faa.gz
gapseq find -p all -b 200 -m auto -t auto $modelB.faa.gz

# (2) Transporter prediction
gapseq find-transport -b 200 $modelA.faa.gz 
gapseq find-transport -b 200 $modelB.faa.gz

# (3) Building Draft Model - based on Reaction-, Pathway-, and Transporter prediction
gapseq draft -r $modelA-all-Reactions.tbl -t $modelA-Transporter.tbl -p $modelA-all-Pathways.tbl -u 200 -l 100 -c $modelA.faa.gz
gapseq draft -r $modelB-all-Reactions.tbl -t $modelB-Transporter.tbl -p $modelB-all-Pathways.tbl -u 200 -l 100 -c $modelB.faa.gz

# (4) Gapfilling
gapseq fill -m $modelA-draft.RDS -n milk.csv -c $modelA-rxnWeights.RDS -g $modelA-rxnXgenes.RDS -b 100
gapseq fill -m $modelB-draft.RDS -n milk.csv -c $modelB-rxnWeights.RDS -g $modelB-rxnXgenes.RDS -b 100
```



##### FBA and FVA prediction of metabolic by-products

Here, we will use the R-Package `sybil` ([Gelius-Dietrich *et al.* (2013) BMC Syst Biol](https://doi.org/10.1186/1752-0509-7-125)) to perform Flux-Balance-Analysis (FBA) and Flux-Variability-Analysis (FVA) with the two reconstructed network models.

First, we define a function, that automatically performs FBA with the minimalization of total flux (MTF) as secondary objective and FVA for all exchange reactions. The function also summarizes the results in a sorted `data.table`.

```R
getMetaboliteProduction <- function(mod) {
  require(sybil)
  require(data.table)
  
  # MTF
  sol.mtf <- optimizeProb(mod, algorithm = "mtf")
  dt.mtf  <- data.table(ex = mod@react_id,
                        mtf.flux = sol.mtf@fluxdist@fluxes[1:mod@react_num])
  dt.mtf.tmp <- copy(dt.mtf[grepl("^EX_cpd[0-9]+_e0", ex)])
  
  # FVA
  sol.fv <- fluxVar(mod, react = mod@react_id[grep("^EX_cpd[0-9]+_e0", mod@react_id)])
  
  dt <- data.table(ex       = rep(mod@react_id[grep("^EX_cpd[0-9]+_e0", mod@react_id)],2),
                   rxn.name = rep(mod@react_name[grep("^EX_cpd[0-9]+_e0", mod@react_id)],2),
                   dir      = c(rep("l",length(grep("^EX_cpd[0-9]+_e0", mod@react_id))),
                                rep("u",length(grep("^EX_cpd[0-9]+_e0", mod@react_id)))),
                   fv       = sol.fv@lp_obj)
  dt <- dcast(dt, ex + rxn.name ~ dir, value.var = "fv")[(u>1e-6 & l >= 0)]
  
  dt <- merge(dt, dt.mtf, by = "ex")
  
  return(dt[order(-mtf.flux)])
}
```

Now, we can apply this function to the network models of *L. delbrueckii* and *S. thermophilus* to predict the top 10 produced metabolic by-products.

```R
library(sybil)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI") # (optional)

ld <- readRDS("ldel.RDS") # for L. delbrueckii
st <- readRDS("sthe.RDS") # for S. thermophilus

getMetaboliteProduction(ld)[1:10]
getMetaboliteProduction(st)[1:10]
```

Output for *L. delbrueckii* (`ld`):

```
                ex               rxn.name           l           u    mtf.flux
 1: EX_cpd00221_e0  D-Lactate-e0 Exchange 0.000000000 4.700648111 4.630126137
 2: EX_cpd00067_e0         H+-e0 Exchange 4.489895516 4.581828957 4.546847983
 3: EX_cpd00108_e0  Galactose-e0 Exchange 2.499994262 2.500000000 2.500000000
 4: EX_cpd00011_e0        CO2-e0 Exchange 0.245315421 0.337271815 0.280318515
 5: EX_cpd00239_e0        H2S-e0 Exchange 0.093096142 0.093119095 0.093119080
 6: EX_cpd00130_e0   L-Malate-e0 Exchange 0.017063114 0.044365180 0.026857932
 7: EX_cpd00029_e0    Acetate-e0 Exchange 0.004859881 0.035769093 0.022378268
 8: EX_cpd00141_e0 Propionate-e0 Exchange 0.007424155 0.007443787 0.007441387
 9: EX_cpd00036_e0  Succinate-e0 Exchange 0.005246866 0.005287035 0.005269831
10: EX_cpd00047_e0    Formate-e0 Exchange 0.000000000 0.003589257 0.003586849
```

Output for S. thermophilus (`st`):

```
                ex                          rxn.name           l            u     mtf.flux
 1: EX_cpd00221_e0             D-Lactate-e0 Exchange 0.000000000 9.405942e+00 7.2350156625
 2: EX_cpd00011_e0                   CO2-e0 Exchange 0.000000000 1.001951e+01 0.4805519366
 3: EX_cpd00029_e0               Acetate-e0 Exchange 0.000000000 8.060152e-01 0.1185652022
 4: EX_cpd00239_e0                   H2S-e0 Exchange 0.088819140 9.015105e-02 0.0888196466
 5: EX_cpd00036_e0             Succinate-e0 Exchange 0.011037101 8.630227e-01 0.0476673792
 6: EX_cpd00363_e0               Ethanol-e0 Exchange 0.000000000 9.405942e+00 0.0019970957
 7: EX_cpd01981_e0 5-Methylthio-D-ribose-e0 Exchange 0.000665698 6.656986e-04 0.0006656986
 8: EX_cpd00020_e0              Pyruvate-e0 Exchange 0.000000000 1.209023e+00 0.0000000000
 9: EX_cpd00071_e0          Acetaldehyde-e0 Exchange 0.000000000 1.209023e+00 0.0000000000
10: EX_cpd00100_e0              Glycerol-e0 Exchange 0.000000000 3.359532e-03 0.0000000000
```

As expected, both organisms produce Lactate in the FBA+MTF solution. The FVA further predicted a lower bound for L-Lactate or D-Lactate production of zero. This is due to the fact, that the models harbour also the capability to produce the respective other Lactate enantiomer. In contrast to *S. thermophilus*, the FBA simulation predicted a release of Galactose by *L. debrueckii*. In fact, *L. debrueckii* is usually reported to be Galactose negative; i.e. does not produce acid from this hexose ([https://bacdive.dsmz.de/strain/6449](https://bacdive.dsmz.de/strain/6449)) and utilized only the glucose part of lactose, while *S. thermophilus* has been reported to be Galactose positive ([https://bacdive.dsmz.de/strain/14786](https://bacdive.dsmz.de/strain/14786)).
