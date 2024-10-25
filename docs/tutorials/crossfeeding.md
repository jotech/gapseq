# Cross-feeding of two gastrointestinal bacteria

> <mark>Please note</mark>: This tutorial assumes a gapseq  <= v1.3.1. An updated version of the tutorial for gapseq >= v1.4.0 is under construction.

##### Background

The intestinal bacterium *Eubacterium rectale* is known to be able to use acetate as energy source under anaerobic conditions and thereby forms butyrate as end product ([Rivère *et al.* (2015) Appl Envrion Microbiol](https://pubmed.ncbi.nlm.nih.gov/26319874/)). Acetate is a common fermentation end product in a number of different other intestinal bacteria, including Bifidobacteria (e.g. *Bifidobacterium longum*). In this tutorial, genome-scale models for *E. rectale* and *B. longum* are reconstructed using **gapseq**. Subsequently, the two models are simulated in co-growth and their interaction is investigated.

*NOTE: All intermediate files produced by the commands below are stored at the github repository (https://github.com/Waschina/gapseq.tutorial.data), which you could download/clone if you wish to start not at the beginning but at a later step of this tutorial.*

##### Input

- Genome assemblies:

  - *Eubacterium rectale* ATCC 33656

    RefSeq: `GCF_000020605.1`

  - *Bifidobacterium longum* NCC2705: 

    RefSeq: `GCF_000007525.1`

- Growth media file: `gf_medium.csv` 

  This is basically a glucose and acetate minimal medium. No amino acids are added since both organisms are likely prototrophic for all proteinogenic amino acids, based on predictions using [GapMind](http://papers.genomics.lbl.gov/cgi-bin/gapView.cgi) ([Price et al. (2019) mSystems](https://doi.org/10.1101/741918 )).  

  E. rectale: [View Gapmind results](http://papers.genomics.lbl.gov/cgi-bin/gapView.cgi?orgs=NCBI__GCF_000020605.1&set=aa); B. longum: [View Gapmind results](http://papers.genomics.lbl.gov/cgi-bin/gapView.cgi?orgs=NCBI__GCF_000007525.1&set=aa)



##### Preparations

Download genome assemblies and gapfill medium. Renaming files.

```sh
#!/bin/bash

# Download genome assemblies 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/605/GCF_000020605.1_ASM2060v1/GCF_000020605.1_ASM2060v1_protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/525/GCF_000007525.1_ASM752v1/GCF_000007525.1_ASM752v1_protein.faa.gz

# Download gapfill-medium file
wget https://github.com/Waschina/gapseq.tutorial.data/raw/master/CF_eure_bilo/gf_medium.csv

# Rename genomes to "eure" (E. rectale) and "bilo" (B. longum) 
mv GCF_000020605.1_ASM2060v1_protein.faa.gz eure.faa.gz
mv GCF_000007525.1_ASM752v1_protein.faa.gz bilo.faa.gz
```



##### gapseq reconstruction 

Now we have the genome sequences and a gapfill medium. That is all we need. Lets reconstruct models:

```sh
#!/bin/bash

modelA="eure"
modelB="bilo"

# Reaction & Pathway prediction
gapseq find -p all -b 200 -m Bacteria $modelA.faa.gz
gapseq find -p all -b 200 -m Bacteria $modelB.faa.gz

# Transporter prediction
gapseq find-transport -b 200 $modelA.faa.gz 
gapseq find-transport -b 200 $modelB.faa.gz

# Building Draft Model - based on Reaction-, Pathway-, and Transporter prediction
gapseq draft -r $modelA-all-Reactions.tbl -t $modelA-Transporter.tbl -p $modelA-all-Pathways.tbl -c $modelA.faa.gz -u 200 -l 100
gapseq draft -r $modelB-all-Reactions.tbl -t $modelB-Transporter.tbl -p $modelB-all-Pathways.tbl -c $modelB.faa.gz -u 200 -l 100

# Gapfilling
gapseq fill -m $modelA-draft.RDS -n gf_medium.csv -c $modelA-rxnWeights.RDS -g $modelA-rxnXgenes.RDS -b 100
gapseq fill -m $modelB-draft.RDS -n gf_medium.csv -c $modelB-rxnWeights.RDS -g $modelB-rxnXgenes.RDS -b 100
```

The final models are stored as R-Object files: `eure.RDS` and `bilo.RDS`, which can be loaded in R using the `readRDS()` command. 


##### Community simulation

Here, we will use the R-Package `BacArena` to perform an agent-based simulation for the co-growth of *B. longum* and *E. rectale*. The following code-block shows the R-source code for a simple community metabolism simulation.

```R
# Load R-packages
library(BacArena)
library(data.table)

# Load reconstructed models
er <- readRDS("eure.RDS") # E. rectale
bl <- readRDS("bilo.RDS") # B. longum

# Small fix to D/L-Lactate secretion (*) and model names
bl <- rmReact(bl, react = "EX_cpd00221_e0")
er@mod_desc <- "E. rectale"
bl@mod_desc <- "B. longum"

# Construct the organism objects for BacArena simulations
eure <- Bac(er)
bilo <- Bac(bl)

# Construct the arena size 10x10 grid cells
arena <- Arena(n = 10, m = 10)

# For each organism, populate randomly 2 grid cells in the Arena as 
# 'starter culture'
arena <- addOrg(arena, eure, amount = 2)
arena <- addOrg(arena, bilo, amount = 2)

# add substrates to arena
arena_subs <- fread("gf_medium.csv") # same as gapfill medium
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]

arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
# Remove acetate from initial substrate list to see effect of Cross-Feeding
arena <- rmSubs(arena, mediac = "EX_cpd00029_e0") 

```

( * *gapseq frequently predicts, that if the organism is producing lactate as fermentation end product,  the optimal solution could involve both  enantiomers: D-/L-Lactate. For plotting & analysis reasons we prohibit the production of D-Lactate to ensure that we see the level of produced Lactate only as one metabolite: L-Lactate. This has otherwise no effect on simulation results.*)

Now we are ready to perform the actual community simulation and plot results:

```R
# Simulation for 13 time steps
CF_sim <- simEnv(arena,time=13, sec_obj = "mtf")

# Plot levels of Acetate, Buyrate, and Lactate as well as growth
par(mfrow=c(1,2))
plotCurves2(CF_sim,legendpos = "topleft",
            subs = c("cpd00211_e0","cpd00029_e0","cpd00159_e0"),
            dict = list(cpd00211_e0 = "Butyrate", 
                        cpd00029_e0 = "Acetate", 
                        cpd00159_e0 = "Lactate"))
```

![](https://github.com/Waschina/gapseq.tutorial.data/raw/master/CF_eure_bilo/CF_eure_bilo.png)

The simulations predicted, that acetate, butyrate, and lactate are produced during co-growth of *E. rectale* and *B. longum*.

Next, let's see, if some of the fermentation products are partially consumed by one of the organisms. This involves a few lines of data wrangling:

```R
# Lets get the exchange fluxs for each grid cell at time step 11
dt.cf <- CF_sim@exchangeslist[[11]]

# Limit output to Acetate, Butyrate, and Lactate
dt.cf <- as.data.table(dt.cf[,c("species",
                                "EX_cpd00029_e0",
                                "EX_cpd00211_e0",
                                "EX_cpd00159_e0")])

# Rename column names (this is just aestetics)
dt.cf <- dt.cf[,.(species, 
                  Acetate = EX_cpd00029_e0, 
                  Butyrate = EX_cpd00211_e0,
                  Lactate = EX_cpd00159_e0)]

# Wide-To-Long table transformation
dt.cf <- melt(dt.cf, 
              id.vars = "species", 
              variable.name = "Metabolite", 
              value.name = "Flux")
dt.cf <- dt.cf[!is.na(Flux)] # rm NA entries (no exchange reaction in grid cell)

# Sum exchanges for each species and metabolite over all 100 grid cells
dt.cf[, .(summed.Flux = sum(Flux)), by = c("species","Metabolite")]
```

The output:

```
      species Metabolite summed.Flux
1: E. rectale    Acetate  -1066.7647
2:  B. longum    Acetate   1925.5036
3: E. rectale   Butyrate   2694.5406
4:  B. longum   Butyrate      0.0000
5:  B. longum    Lactate    806.4733
```

We can see, that *B. longum* (bilo) secretes acetate and lactate as main end product. Approximately 55 % of acetate is consumed by *E. rectale* (eure). *B. longum* produces in addition lactate and *E. rectale* secretes butyrate.

