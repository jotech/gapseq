# Lactic acid bacteria in yogurt production

*Lactobacillus delbrueckii* subsp. *bulgaricus* and *Streptococcus thermophilus* are both lactic acid bacteria, which are frequently used for the production of yogurt. The organisms differ in their lactose degradation and their fermentation end products. In this tutorial, genome-scale metabolic models will be reconstructed using **gapseq**.



##### Input files

- Genome assemblies:
  - *Lactobacillus delbrueckii* subsp. *bulgaricus* ATCC 11842 = JCM 1002. RefSeq: `GCF_000056065.1`
  - *Streptococcus thermophilus* ATCC 19258. RefSeq: `GCF_010120595.1`
- Growth media (milk) file: `milk.csv` 
  This growth media is based on the main ingredients described for whole milk (without fatty acids), as listed [here](https://frida.fooddata.dk/food/1265?lang=en).

*NOTE: All intermediate files produced by the commands below are stored at the github repository [gapseq.tutorial.data](https://github.com/Waschina/gapseq.tutorial.data), which you could  download/clone if you wish to start not at the beginning but at a later  step of this tutorial.*



##### Preparations

Download required files and rename assembly filed for the ease of this tutorial.

```sh
mkdir yoghurt
cd yoghurt

# Download genome assemblies from NCBI (RefSeq)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/056/065/GCF_000056065.1_ASM5606v1/GCF_000056065.1_ASM5606v1_genomic.fna.gz .
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/010/120/595/GCF_010120595.1_ASM1012059v1/GCF_010120595.1_ASM1012059v1_genomic.fna.gz .

# Download gapfill-medium file from github
wget https://github.com/Waschina/gapseq.tutorial.data/raw/master/yogurt/milk.csv .

# Rename genomes
mv GCF_000056065.1_ASM5606v1_genomic.fna.gz ldel.fna.gz
mv GCF_010120595.1_ASM1012059v1_genomic.fna.gz sthe.fna.gz
```



##### gapseq reconstruction pipeline

(1) Reaction & pathway prediction
(2) Transporter prediction
(3) Draft model reconstruction
(4) Gapfilling

```sh
modelA="ldel"
modelB="sthe"

# If not set already (e.g via .bashrc): set the path to gapseq
# There are different ways to do this. One example:
gapseq=~/path/to/gapseq/./gapseq

# (1) Reaction & Pathway prediction
$gapseq find -p all -b 200 $modelA.fna.gz
$gapseq find -p all -b 200 $modelB.fna.gz

# (2) Transporter prediction
$gapseq find-transport -b 200 $modelA.fna.gz 
$gapseq find-transport -b 200 $modelB.fna.gz

# (3) Building Draft Model - based on Reaction-, Pathway-, and Transporter prediction
$gapseq draft -r $modelA-all-Reactions.tbl -t $modelA-Transporter.tbl -p $modelA-all-Pathways.tbl -c $modelA.fna.gz -u 200 -l 100
$gapseq draft -r $modelB-all-Reactions.tbl -t $modelB-Transporter.tbl -p $modelB-all-Pathways.tbl -c $modelB.fna.gz -u 200 -l 100

# (4) Gapfilling
$gapseq fill -m $modelA-draft.RDS -n milk.csv -c $modelA-rxnWeights.RDS -g $modelA-rxnXgenes.RDS -b 100
$gapseq fill -m $modelB-draft.RDS -n milk.csv -c $modelB-rxnWeights.RDS -g $modelB-rxnXgenes.RDS -b 100
```



##### FBA and FVA prediction of metabolic by products

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
ld <- readRDS("ldel.RDS") # for L. delbrueckii
st <- readRDS("sthe.RDS") # for S. thermophilus

getMetaboliteProduction(ld)[1:10]
getMetaboliteProduction(st)[1:10]
```

Output for *L. delbrueckii* (`ld`):

```
                ex                      rxn.name           l            u    mtf.flux
 1: EX_cpd00067_e0                H+-e0 Exchange 5.251791165 10.269349014 5.306939206
 2: EX_cpd00159_e0         L-Lactate-e0 Exchange 0.000000000  9.725953947 5.277189818
 3: EX_cpd00108_e0         Galactose-e0 Exchange 2.866724500  5.000000000 2.917862002
 4: EX_cpd00011_e0               CO2-e0 Exchange 0.000000000  1.775317874 0.255635652
 5: EX_cpd00036_e0         Succinate-e0 Exchange 0.000000000  3.643758830 0.094979517
 6: EX_cpd00141_e0        Propionate-e0 Exchange 0.064777627  0.065703755 0.064778207
 7: EX_cpd00029_e0           Acetate-e0 Exchange 0.000000000  3.286683325 0.018148342
 8: EX_cpd00281_e0              GABA-e0 Exchange 0.000000000  0.078770095 0.006628222
 9: EX_cpd03091_e0 5'-Deoxyadenosine-e0 Exchange 0.002416292  0.002416363 0.002416301
10: EX_cpd00239_e0               H2S-e0 Exchange 0.001610510  0.092611584 0.001610867
```

Output for S. thermophilus (`st`):

```
                ex                      rxn.name            l            u    mtf.flux
 1: EX_cpd00159_e0         L-Lactate-e0 Exchange 0.0000000000 19.373109318 6.145385831
 2: EX_cpd00047_e0           Formate-e0 Exchange 0.0000000000 16.347297990 1.433212770
 3: EX_cpd00029_e0           Acetate-e0 Exchange 0.0000000000 23.775767754 1.108637403
 4: EX_cpd00036_e0         Succinate-e0 Exchange 0.0000000000 11.908553089 1.101411232
 5: EX_cpd01947_e0              BDOH-e0 Exchange 0.0000000000 10.586611615 0.516915783
 6: EX_cpd00239_e0               H2S-e0 Exchange 0.0850060915  0.087118626 0.087118568
 7: EX_cpd00141_e0        Propionate-e0 Exchange 0.0006958976  9.083687997 0.041328170
 8: EX_cpd00363_e0           Ethanol-e0 Exchange 0.0000000000 19.373109318 0.028435080
 9: EX_cpd00100_e0          Glycerol-e0 Exchange 0.0000000000  0.004556206 0.004556206
10: EX_cpd03091_e0 5'-Deoxyadenosine-e0 Exchange 0.0042126479  0.004213093 0.004212667
```

As expected, both organisms produce Lactate in the FBA+MTF solution. The FVA further predicted a lower bound for L-Lactate production of zero. This is due to the fact, that the models harbour also the capability to produce the enantiomer D-Lactate. In contrast to *S. thermophilus*, the FBA simulation predicted a release of Galactose by *L. debrueckii*. In fact, *L. debrueckii* is usually reported to be Galactose negative; i.e. does not produce acid from this hexose ([https://bacdive.dsmz.de/strain/6449](https://bacdive.dsmz.de/strain/6449)), while *S. thermophilus* has been reported to be Galactose positive ([https://bacdive.dsmz.de/strain/14786](https://bacdive.dsmz.de/strain/14786)).
