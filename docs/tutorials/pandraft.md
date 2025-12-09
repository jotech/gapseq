# Example Workflow for pan-draft
2025-12-04

### Background

This example demonstrates how to use **pan-draft** with gapseq to reconstruct a pan-genome-scale metabolic model. We use genomes (mainly MAGs) of the species *Blautia hansenii* from the Unified Human Gastrointestinal Genome (UHGG) version 2.0.2 catalog as an example.

1. Download genome metadata and sequences from UHGG.
2. Extract nucleotide sequences from GFF files.
3. Generate individual draft genome-scale metabolic models with gapseq.
4. Build a pan-draft model using pan-draft.
5. Gapfill the pan-draft model to obtain a functional model as representation of the *Blautia hansenii*-representative metabolism.

#### 1a. Download UHGG Genome Metadata

We begin by retrieving the UHGG catalog, which contains metadata and download links for all genomes.


``` sh
# get UHGG catalog (v2.0.2)
wget -q https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/genomes-all_metadata.tsv
```

#### 1b. Identify Available Genomes for *Blautia hansenii*

Next, we check how many genomes in the catalog belong to *Blautia hansenii*.


``` sh
# Check how many genomes are in the catalog for Blautia hansenii
awk -F'\t' '$15 ~ /Blautia hansenii/ {count++} END {print count}' genomes-all_metadata.tsv
```

```
## 117
```
This returns the number of genomes available for the target species. (117)

#### 1c. Download Genome Files

We now download the GFF files for all *Blautia hansenii* genomes.
These GFF files include the annotations, and importantly, the genomic nucleotide sequences at the end of the file (after the `##FASTA` line).


``` sh
# Download genomes for species Blautia hansenii
mkdir -p gff
awk -F'\t' '$15 ~ /Blautia hansenii/ {print $20}' genomes-all_metadata.tsv \
  | while read url; do
      wget -nc -P gff "$url"
    done
```

#### 2. Extract Nucleotide Sequences

The nucleotide sequences are embedded within the GFF files.
We extract them and store them in separate FASTA files, which will serve as the input for gapseq.


``` sh
# extract genomic nucleotide sequences and save them as fasta
mkdir -p fna
for gff in gff/*.gff.gz; do
    base=$(basename "$gff" .gff.gz)
    zcat "$gff" | awk '/^##FASTA$/ {found=1; next} found' > "fna/${base}.fna"
done
```

#### 3. Reconstruct Draft Metabolic Networks with gapseq

For each genome, we now reconstruct an individual draft metabolic model using gapseq.
This includes:

- Identifying metabolic genes (`gapseq find`).
- Detecting transporters (`gapseq find-transport`).
- Drafting the metabolic network (`gapseq draft`).

<mark>Please note:</mark> This step takes quite some time to complete. To optimize the process, we use [GNU parallel](https://www.gnu.org/software/parallel/) in the following code, to reconstruct four draft models in parallel.


``` sh
# the fasta files are now the input for gapseq:
mkdir -p models

parallel -j 4 --bar '
    fna={}
    id=$(basename "$fna" .fna)
    echo "Reconstructing draft model for $id"

    gapseq find -p all -t Bacteria -m Bacteria -f models -b 200 -A diamond "$fna"
    gapseq find-transport -f models -b 200 -A diamond "$fna"
    gapseq draft -f models -u 200 -l 100 \
        -r models/${id}-all-Reactions.tbl \
        -t models/${id}-Transporter.tbl \
        -p models/${id}-all-Pathways.tbl
' ::: fna/*.fna
```

At this point, we have a set of draft metabolic models, one per genome.

#### 4. Build the Pan-Draft Model

Now, we combine the individual draft metabolic models into a pan-draft model.
This model captures the core metabolic reactions shared across genomes, as well as accessory reactions.


``` sh
# reconstruct pan model
mkdir -p pan_model
gapseq pan -m models/ -w models/ -f pan_model/
```

```
## 
## Loading input data, be patient... 
## Reading input from directory:	 models/ 
## models//MGYG000001704-draft.RDS models//MGYG000004911-draft.RDS models//MGYG000009596-draft.RDS models//MGYG000014854-draft.RDS models//MGYG000018898-draft.RDS models//MGYG000018936-draft.RDS models//MGYG000019565-draft.RDS models//MGYG000026226-draft.RDS models//MGYG000027117-draft.RDS models//MGYG000029914-draft.RDS models//MGYG000033847-draft.RDS models//MGYG000037111-draft.RDS models//MGYG000043121-draft.RDS models//MGYG000043600-draft.RDS models//MGYG000052579-draft.RDS models//MGYG000054967-draft.RDS models//MGYG000055367-draft.RDS models//MGYG000056559-draft.RDS models//MGYG000062338-draft.RDS models//MGYG000062734-draft.RDS models//MGYG000064419-draft.RDS models//MGYG000064946-draft.RDS models//MGYG000067028-draft.RDS models//MGYG000069312-draft.RDS models//MGYG000069581-draft.RDS models//MGYG000069602-draft.RDS models//MGYG000069775-draft.RDS models//MGYG000071151-draft.RDS models//MGYG000071376-draft.RDS models//MGYG000076053-draft.RDS models//MGYG000076749-draft.RDS models//MGYG000078674-draft.RDS models//MGYG000082887-draft.RDS models//MGYG000083382-draft.RDS models//MGYG000084533-draft.RDS models//MGYG000087186-draft.RDS models//MGYG000088773-draft.RDS models//MGYG000090503-draft.RDS models//MGYG000099188-draft.RDS models//MGYG000107300-draft.RDS models//MGYG000108542-draft.RDS models//MGYG000115300-draft.RDS models//MGYG000116616-draft.RDS models//MGYG000119248-draft.RDS models//MGYG000122398-draft.RDS models//MGYG000123443-draft.RDS models//MGYG000123847-draft.RDS models//MGYG000131877-draft.RDS models//MGYG000142364-draft.RDS models//MGYG000143354-draft.RDS models//MGYG000144227-draft.RDS models//MGYG000145050-draft.RDS models//MGYG000145518-draft.RDS models//MGYG000145816-draft.RDS models//MGYG000147852-draft.RDS models//MGYG000155173-draft.RDS models//MGYG000155379-draft.RDS models//MGYG000159398-draft.RDS models//MGYG000160973-draft.RDS models//MGYG000163223-draft.RDS models//MGYG000168116-draft.RDS models//MGYG000177827-draft.RDS models//MGYG000180265-draft.RDS models//MGYG000182417-draft.RDS models//MGYG000184690-draft.RDS models//MGYG000185155-draft.RDS models//MGYG000188105-draft.RDS models//MGYG000188348-draft.RDS models//MGYG000188358-draft.RDS models//MGYG000189392-draft.RDS models//MGYG000190963-draft.RDS models//MGYG000193789-draft.RDS models//MGYG000194008-draft.RDS models//MGYG000194762-draft.RDS models//MGYG000199243-draft.RDS models//MGYG000204172-draft.RDS models//MGYG000205451-draft.RDS models//MGYG000206027-draft.RDS models//MGYG000211179-draft.RDS models//MGYG000214161-draft.RDS models//MGYG000214917-draft.RDS models//MGYG000215554-draft.RDS models//MGYG000215948-draft.RDS models//MGYG000219677-draft.RDS models//MGYG000227318-draft.RDS models//MGYG000230891-draft.RDS models//MGYG000231352-draft.RDS models//MGYG000232767-draft.RDS models//MGYG000234137-draft.RDS models//MGYG000236184-draft.RDS models//MGYG000236459-draft.RDS models//MGYG000236698-draft.RDS models//MGYG000238827-draft.RDS models//MGYG000239425-draft.RDS models//MGYG000240689-draft.RDS models//MGYG000248090-draft.RDS models//MGYG000248131-draft.RDS models//MGYG000249630-draft.RDS models//MGYG000250009-draft.RDS models//MGYG000261262-draft.RDS models//MGYG000262284-draft.RDS models//MGYG000262314-draft.RDS models//MGYG000262907-draft.RDS models//MGYG000270086-draft.RDS models//MGYG000273714-draft.RDS models//MGYG000276051-draft.RDS models//MGYG000277845-draft.RDS models//MGYG000278797-draft.RDS models//MGYG000280146-draft.RDS models//MGYG000280354-draft.RDS models//MGYG000281077-draft.RDS models//MGYG000281572-draft.RDS models//MGYG000281755-draft.RDS models//MGYG000285258-draft.RDS models//MGYG000285443-draft.RDS models//MGYG000288642-draft.RDS models//MGYG000288705-draft.RDS 
## Reading input from directory:	 models/ 
## models//MGYG000001704-all-Pathways.tbl models//MGYG000004911-all-Pathways.tbl models//MGYG000009596-all-Pathways.tbl models//MGYG000014854-all-Pathways.tbl models//MGYG000018898-all-Pathways.tbl models//MGYG000018936-all-Pathways.tbl models//MGYG000019565-all-Pathways.tbl models//MGYG000026226-all-Pathways.tbl models//MGYG000027117-all-Pathways.tbl models//MGYG000029914-all-Pathways.tbl models//MGYG000033847-all-Pathways.tbl models//MGYG000037111-all-Pathways.tbl models//MGYG000043121-all-Pathways.tbl models//MGYG000043600-all-Pathways.tbl models//MGYG000052579-all-Pathways.tbl models//MGYG000054967-all-Pathways.tbl models//MGYG000055367-all-Pathways.tbl models//MGYG000056559-all-Pathways.tbl models//MGYG000062338-all-Pathways.tbl models//MGYG000062734-all-Pathways.tbl models//MGYG000064419-all-Pathways.tbl models//MGYG000064946-all-Pathways.tbl models//MGYG000067028-all-Pathways.tbl models//MGYG000069312-all-Pathways.tbl models//MGYG000069581-all-Pathways.tbl models//MGYG000069602-all-Pathways.tbl models//MGYG000069775-all-Pathways.tbl models//MGYG000071151-all-Pathways.tbl models//MGYG000071376-all-Pathways.tbl models//MGYG000076053-all-Pathways.tbl models//MGYG000076749-all-Pathways.tbl models//MGYG000078674-all-Pathways.tbl models//MGYG000082887-all-Pathways.tbl models//MGYG000083382-all-Pathways.tbl models//MGYG000084533-all-Pathways.tbl models//MGYG000087186-all-Pathways.tbl models//MGYG000088773-all-Pathways.tbl models//MGYG000090503-all-Pathways.tbl models//MGYG000099188-all-Pathways.tbl models//MGYG000107300-all-Pathways.tbl models//MGYG000108542-all-Pathways.tbl models//MGYG000115300-all-Pathways.tbl models//MGYG000116616-all-Pathways.tbl models//MGYG000119248-all-Pathways.tbl models//MGYG000122398-all-Pathways.tbl models//MGYG000123443-all-Pathways.tbl models//MGYG000123847-all-Pathways.tbl models//MGYG000131877-all-Pathways.tbl models//MGYG000142364-all-Pathways.tbl models//MGYG000143354-all-Pathways.tbl models//MGYG000144227-all-Pathways.tbl models//MGYG000145050-all-Pathways.tbl models//MGYG000145518-all-Pathways.tbl models//MGYG000145816-all-Pathways.tbl models//MGYG000147852-all-Pathways.tbl models//MGYG000155173-all-Pathways.tbl models//MGYG000155379-all-Pathways.tbl models//MGYG000159398-all-Pathways.tbl models//MGYG000160973-all-Pathways.tbl models//MGYG000163223-all-Pathways.tbl models//MGYG000168116-all-Pathways.tbl models//MGYG000177827-all-Pathways.tbl models//MGYG000180265-all-Pathways.tbl models//MGYG000182417-all-Pathways.tbl models//MGYG000184690-all-Pathways.tbl models//MGYG000185155-all-Pathways.tbl models//MGYG000188105-all-Pathways.tbl models//MGYG000188348-all-Pathways.tbl models//MGYG000188358-all-Pathways.tbl models//MGYG000189392-all-Pathways.tbl models//MGYG000190963-all-Pathways.tbl models//MGYG000193789-all-Pathways.tbl models//MGYG000194008-all-Pathways.tbl models//MGYG000194762-all-Pathways.tbl models//MGYG000199243-all-Pathways.tbl models//MGYG000204172-all-Pathways.tbl models//MGYG000205451-all-Pathways.tbl models//MGYG000206027-all-Pathways.tbl models//MGYG000211179-all-Pathways.tbl models//MGYG000214161-all-Pathways.tbl models//MGYG000214917-all-Pathways.tbl models//MGYG000215554-all-Pathways.tbl models//MGYG000215948-all-Pathways.tbl models//MGYG000219677-all-Pathways.tbl models//MGYG000227318-all-Pathways.tbl models//MGYG000230891-all-Pathways.tbl models//MGYG000231352-all-Pathways.tbl models//MGYG000232767-all-Pathways.tbl models//MGYG000234137-all-Pathways.tbl models//MGYG000236184-all-Pathways.tbl models//MGYG000236459-all-Pathways.tbl models//MGYG000236698-all-Pathways.tbl models//MGYG000238827-all-Pathways.tbl models//MGYG000239425-all-Pathways.tbl models//MGYG000240689-all-Pathways.tbl models//MGYG000248090-all-Pathways.tbl models//MGYG000248131-all-Pathways.tbl models//MGYG000249630-all-Pathways.tbl models//MGYG000250009-all-Pathways.tbl models//MGYG000261262-all-Pathways.tbl models//MGYG000262284-all-Pathways.tbl models//MGYG000262314-all-Pathways.tbl models//MGYG000262907-all-Pathways.tbl models//MGYG000270086-all-Pathways.tbl models//MGYG000273714-all-Pathways.tbl models//MGYG000276051-all-Pathways.tbl models//MGYG000277845-all-Pathways.tbl models//MGYG000278797-all-Pathways.tbl models//MGYG000280146-all-Pathways.tbl models//MGYG000280354-all-Pathways.tbl models//MGYG000281077-all-Pathways.tbl models//MGYG000281572-all-Pathways.tbl models//MGYG000281755-all-Pathways.tbl models//MGYG000285258-all-Pathways.tbl models//MGYG000285443-all-Pathways.tbl models//MGYG000288642-all-Pathways.tbl models//MGYG000288705-all-Pathways.tbl 
## loading completed 
## 
## The sizes of input lists are consistent, the number of loaded model is 117 
## 
## 
## The total # of rxn is: 2026 
## The # of strict core rxn (all mod) is: 205 
## The # of core rxn ( >= 95 % mod) is: 213 
## The # of shell rxn ( >= 5 % mod) is: 1185 
## The # of cloud rxn ( < 5 % mod) is: 423 
## 
## Let's standardize the name of the duplicated compounds:
## modifying 5-Methyltetrahydrofolate-c0 
## modifying Pyridoxal phosphate-c0 
## modifying 2-Demethylmenaquinone 8-c0 
## modifying Calomide-c0 
## modifying phosphatidylethanolamine dioctadecanoyl-c0 
## modifying Dianteisoheptadecanoylphosphatidylethanolamine-c0 
## modifying Diisoheptadecanoylphosphatidylethanolamine-c0 
## 
## th:  0.06
## Constructing draft model... 
## Loading required package: stringr
## 	completed
```

#### 5a. Predict the Gapfill Medium

To make the draft model functional, gapseq first predicts a minimal growth medium that supports biomass production in the pan-model.


``` sh
# predict gapfill medium for pan model
gapseq medium -m pan_model/panModel-draft.RDS -p pan_model/panModel-tmp-Pathways.tbl -f pan_model/
```

#### 5b. Gapfill the Pan-Model

Finally, we perform gapfilling, which adds reactions to ensure the model can simulate growth.


``` sh
# gapfill pan model
gapseq fill -m pan_model/panModel-draft.RDS -n pan_model/panModel-medium.csv -b 100
```

### Summary

The central output files are `pan_model.RDS` and `pan_model.xml`, both contain the species-level representative metabolic network of *Blautia hansenii*. The RDS-file is a R-object, that can be loaded in R using `mod <- readRDS("pan_model.RDS")`. The XML-file is a SBML representation of the network model and can be used in most constraint-based modelling software tools.

This workflow illustrates how pan-draft can improve the reconstruction of genome-scale metabolic models for species represented by multiple metagenome-assembled genomes (MAGs), yielding more complete and robust models compared to single-genome approaches. The workflow is not limited by the number of genomes that serve as input.


