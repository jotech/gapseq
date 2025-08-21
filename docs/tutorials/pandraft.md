# Example Workflow for pan-draft

### Background

This example demonstrates how to use **pan-draft** with gapseq to reconstruct a pan-genome-scale metabolic model. We use genomes (mainly MAGs) of the species *Blautia hansenii* from the Unified Human Gastrointestinal Genome (UHGG) version 2.0.2 catalog as an example.

1. Download genome metadata and sequences from UHGG.
2. Extract nucleotide sequences from GFF files.
3. Generate individual draft genome-scale metabolic models with gapseq.
4. Build a pan-draft model using pan-draft.
5. Gapfill the pan-draft model to obtain a functional model as representation of the *Blautia hansenii*-representative metabolism.

#### 1a. Download UHGG Genome Metadata

We begin by retrieving the UHGG catalog, which contains metadata and download links for all genomes.

```sh
# get UHGG catalog (v2.0.2)
wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/genomes-all_metadata.tsv
```

#### 1b. Identify Available Genomes for *Blautia hansenii*

Next, we check how many genomes in the catalog belong to *Blautia hansenii*.

```sh
# Check how many genomes are in the catalog for Blautia hansenii
awk -F'\t' '$15 ~ /Blautia hansenii/ {count++} END {print count}' genomes-all_metadata.tsv
```
This will return the number of genomes available for the target species. (117)

#### 1c. Download Genome Files

We now download the GFF files for all *Blautia hansenii* genomes.
These GFF files include the annotations, and importantly, the genomic nucleotide sequences at the end of the file (after the `##FASTA` line).

```sh
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

```sh
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

<mark>Please note:</mark> This step takes quite some time to complete. However, one could modify the code below to improve computation time through parallelization, e.g., via [GNU parallel](https://www.gnu.org/software/parallel/).

```sh
# the fasta files are now the input for gapseq:
mkdir -p models
for fna in fna/*.fna; do
    id=$(basename $fna .fna)
    echo "Reconstructing draft model for $id"
    gapseq find -p all -t Bacteria -m Bacteria -f models -b 200 -A diamond $fna
    gapseq find-transport -f models -b 200 -A diamond $fna
    gapseq draft -f models -u 200 -l 100 -r models/${id}-all-Reactions.tbl -t models/${id}-Transporter.tbl -p models/${id}-all-Pathways.tbl
done
```

At this point, we have a set of draft metabolic models, one per genome.

#### 4. Build the Pan-Draft Model

Now, we combine the individual draft metabolic models into a pan-draft model.
This model captures the core metabolic reactions shared across genomes, as well as accessory reactions.

```sh
# reconstruct pan model
mkdir -p pan_model
gapseq pan -m models/ -w models/ -f pan_model/
```

#### 5a. Predict the Gapfill Medium

To make the draft model functional, gapseq first predicts a minimal growth medium that supports biomass production in the pan-model.

```sh
# predict gapfill medium for pan model
gapseq medium -m pan_model/panModel-draft.RDS -p pan_model/panModel-tmp-Pathways.tbl -f pan_model/
```

#### 5b. Gapfill the Pan-Model

Finally, we perform gapfilling, which adds reactions to ensure the model can simulate growth.

```sh
# gapfill pan model
gapseq fill -m pan_model/panModel-draft.RDS -n pan_model/panModel-medium.csv
```

### Summary

This workflow illustrates how pan-draft can improve the reconstruction of genome-scale metabolic models for species represented by multiple metagenome-assembled genomes (MAGs), yielding more complete and robust models compared to single-genome approaches. The workflow is not limited by the number of genomes that serve as input.
