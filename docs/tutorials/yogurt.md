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
                ex              rxn.name          l          u   mtf.flux
 1: EX_cpd00221_e0 D-Lactate-e0 Exchange 0.00000000 4.66688283 4.57527357
 2: EX_cpd00067_e0        H+-e0 Exchange 4.41375845 4.56576030 4.53889042
 3: EX_cpd00108_e0 Galactose-e0 Exchange 2.49998944 2.50000000 2.50000000
 4: EX_cpd00011_e0       CO2-e0 Exchange 0.32602433 0.47806843 0.35293482
 5: EX_cpd00371_e0  Propanal-e0 Exchange 0.16765414 0.16773865 0.16769626
 6: EX_cpd00047_e0   Formate-e0 Exchange 0.09438168 0.16513265 0.16509025
 7: EX_cpd00013_e0       NH3-e0 Exchange 0.13880177 0.26199257 0.16445716
 8: EX_cpd00239_e0       H2S-e0 Exchange 0.09297787 0.09302013 0.09302010
 9: EX_cpd00324_e0      MTTL-e0 Exchange 0.08651320 0.08655546 0.08655540
10: EX_cpd00029_e0   Acetate-e0 Exchange 0.07586589 0.17284734 0.07654518
```

Output for S. thermophilus (`st`):

```
                ex                          rxn.name            l            u     mtf.flux
 1: EX_cpd00221_e0             D-Lactate-e0 Exchange 0.0000000000 9.2094468818 8.6588710458
 2: EX_cpd00011_e0                   CO2-e0 Exchange 0.2558426738 9.9070817239 0.5225485126
 3: EX_cpd00239_e0                   H2S-e0 Exchange 0.0868496916 0.0884168597 0.0868510077
 4: EX_cpd00281_e0                  GABA-e0 Exchange 0.0644172061 0.0644198887 0.0644184793
 5: EX_cpd00092_e0                Uracil-e0 Exchange 0.0469340783 0.0469354195 0.0469353179
 6: EX_cpd00036_e0             Succinate-e0 Exchange 0.0023487402 0.3596882237 0.0454283729
 7: EX_cpd00029_e0               Acetate-e0 Exchange 0.0000000000 0.2101233068 0.0268075147
 8: EX_cpd00100_e0              Glycerol-e0 Exchange 0.0000000000 0.0039510794 0.0039510794
 9: EX_cpd00363_e0               Ethanol-e0 Exchange 0.0000000000 9.2094468818 0.0023487447
10: EX_cpd01981_e0 5-Methylthio-D-ribose-e0 Exchange 0.0007829134 0.0007829149 0.0007829149
```

As expected, both organisms produce Lactate in the FBA+MTF solution. The FVA further predicted a lower bound for L-Lactate or D-Lactate production of zero. This is due to the fact, that the models harbour also the capability to produce the respective other Lactate enantiomer. In contrast to *S. thermophilus*, the FBA simulation predicted a release of Galactose by *L. debrueckii*. In fact, *L. debrueckii* is usually reported to be Galactose negative; i.e. does not produce acid from this hexose ([https://bacdive.dsmz.de/strain/6449](https://bacdive.dsmz.de/strain/6449)) and utilized only the glucose part of lactose, while *S. thermophilus* has been reported to be Galactose positive ([https://bacdive.dsmz.de/strain/14786](https://bacdive.dsmz.de/strain/14786)).
