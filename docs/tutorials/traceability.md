# Reaction tracing

Keyword: **"Traceability"**

Automatic metabolic network reconstruction is a multi-step process. Reactions and metabolites are added to a metabolic network model in various steps and based on different criteria. ***gapseq* records for each reaction, why and at which stage it was added to the network.** While this tracing of reactions was originally intended for debugging, we believe that this information is also highly relevant for users to understand their output network and for potential manual refinement / curation efforts. 

gapseq stores the reaction tracing information in the reaction's attributes data.frame of the *cobrar::ModelOrg* R-object, which is saved by gapseq in the *model.RDS* file. Here's a short code-snippet to access the reaction tracing information in R:

```R
library(cobrar)

mod <- readRDS("model.RDS") # adjust with the respective model name

# Retrieve data.frame for reaction attriubutes (it's a big table)
mod@react_attr

# Retrieve only reaction IDs, name, ec and the reactions'  
# sources/origin in the reconstruction process
mod@react_attr[,c("seed", "name", "ec", "gs.origin")]

# In case you wish to save this table to a csv-file:
write.csv(mod@react_attr[,c("seed", "name", "ec", "gs.origin")],
          file = "model_reaction_tracing.csv")
```



**Here's what the codes in the `gs.origin` column refer to:**


`0`. Reaction added due to sequence homology to reference proteins

`1`. Gapfilling Step 1 – Enable flux through biomass reaction / growth

`2`. Gapfilling Step 2 – Biosynthesis of biomass components

`3`. Gapfilling Step 3 – Alternative carbon/energy sources

`4`. Gapfilling Step 4 – Potential metabolic products

`5`. Added transport reaction due to present pertinent pathway/reaction. See below for details.

`6`. Biomass reaction

`7`. Exchange reactions

`8`. Diffusion reactions

`9`. Reaction added due to pathway completion 

`10`. Reaction added after using `./gapseq adapt`



Please note: The reaction tracing information is currently not included in the output SBML-files of the models. This will be included in a future update of *gapseq*.

Addition to code `5`: For instance, if an organism harbors the genes for the [degradation of Tryptophan to Indole-3-Propionic Acid (IPA)](https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=PWY-8017), the transport reaction for IPA is automatically added to the model. This automatism is meant to allow the production of metabolites, where the producing-pathways and the cellular export is known from experiments, but where the transport mechanism (i.e. transporter protein) remains unknown. Currently only a few transporters are added based on this rule:

- Indole-3-propionate (cpd27289) transport if [PWY-8017](https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=PWY-8017) is present
- Phenypropanoate (cpd03343) transport if [PWY-8014](https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=PWY-8014) is present
- Phloretate transport if [PWY-8016](https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=PWY-8016) is present
- Butyrate transport if Butyryl-CoA dehydrogenase (Etf complex)) is present (rxn90001)
