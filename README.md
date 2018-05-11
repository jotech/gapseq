# gapseq
Informed curation and gapfilling of metabolic networks


## Example workflow
1) Needed input files
* Genome (fna)
* Draft model (xml)

2) Seed database corrections
* python src/seedcorrect.py dat/myb71.xml

3) Obtaining candidate reactions
* ./gapseq.sh -p core dat/myb71.fna
* ./gapseq.sh -p degradation dat/myb71.fna
* ./transporter.sh dat/myb71.fna

4) Gapfilling
* Rscript gf.suite.R -m ./dat/myb71_corrected.xml -c ./dat/myb71-core-Reactions.lst,dat/myb71-degradation-Reactions.lst,dat/myb71-Transporter.lst -n dat/media/TSBmed.csv
