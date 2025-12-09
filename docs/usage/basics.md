# Basic usage

gapseq has two basic use cases. On the one hand, it can be used to predict **metabolic pathways** and on the other hand, **metabolic networks** can be reconstructed.
Both use cases are linked in the sense that for the reconstruction of metabolic networks, the predicted pathways and the corresponding reactions serve as input.
The following figure summarizes the workflow.


<p align="center">
<img src="https://github.com/jotech/gapseq/raw/master/docs/gfx/flowchart.png" alt="alt text" title="Title"/>
<i>Overview of basic usage in gapseq</i>
</p>

As a quickstart, all predictions can be done by simply using
```
./gapseq doall toy/myb71.faa.gz
```
## Pathway and transporter prediction
Search for metabolic pathways and transporter:
```
./gapseq find -p all toy/myb71.faa.gz
./gapseq find-transport toy/myb71.faa.gz
```

Search for chitin pathways:
```
./gapseq find -p chitin toy/myb71.faa.gz
```

Check for a certain enzyme availability, for example cytochrome c oxidases (oxidase test)
```
./gapseq find -e 1.9.3.1 toy/ecoli.faa.gz
```
Search for enzymes by name (e.g.: bacterial ligning degradation):
```
./gapseq find -r "dye-decolorizing peroxidase" toy/ecoli.faa.gz
```
(quotation marks are needed for reaction names that contain whitespaces)

## Draft network reconstruction and gapfilling
Creation of draft model
```
 ./gapseq draft -r toy/myb71-all-Reactions.tbl -t toy/myb71-Transporter.tbl -p toy/myb71-all-Pathways.tbl
```
Gap filling
```
./gapseq fill -m toy/myb71-draft.RDS -n dat/media/TSBmed.csv
```
