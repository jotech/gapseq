# Output files
In this section, the various final and intermediate output files of gapseq are going to be described.

## gapseq find
### *-Pathways.tbl
This file contains detailed information for all pathways that were considered by the ``gapseq find`` command.
It is a tab-separated text file with the following columns: 

| ID         | Name         | Prediction                   | Completeness                   | VagueReactions                                  | KeyReactions                            | KeyReactionsFound             | ReactionsFound           |
| -          | -        | -                            | -                              | -                                               | -                                       | -                             | -                        |
| Pathway ID | Pathway Name | Inferred presence true/false | Ratio of found/total reactions (%) | Number of reactions without available sequences | Total number of important key reactions | Number of found key reactions | Names of found reactions |

### *-Reactions.tbl
This tab-separated text file contains detailed information about all checked reactions.

| rxn         | name          | ec        | bihit             | qseqid       | pident                          | evalue       | bitscore  | qcovs                      | stitle        | sstart             | send             | pathway            | status                    | pathway.status                 | dbhit                  | complex                  | exception                      | complex.status            |
| -           | -             | -         | -                 | -            | -                               | -            | -         | -                          | -             | -                  | -                | -                  | -                         | -                              | -                      | -                        | -                              | -                         |
| Reaction ID | Reaction name | EC number | Bidirectional hit | Query Seq-id | Percentage of identical matches | Expect value | Bit score | Query Coverage Per Subject | Subject Title | Start of alignment | End of alignment | Associated pathway | Blast status of reaction* | Status of associated pathway** | Mapped model reactions | Detected protein complex | Higher identity cutoff used*** | Status of protein complex |

* *The blast status of a reaction informs about the result of the homology search. It is defined to be: ``bad_blast`` (blast hit with lower quality, i.e. lower bitscore, coverage, or identity than needed for cutoffs), ``good_blast`` (all cutoffs satisfying blast hit), ``no_blast`` (no blast hit found), ``no_seq_data`` (no sequence data available), ``spontaneous`` (no enzyme needed).
* **The status of the associated pathway provides background to the criteria by which a pathway was predicted. The following values are possible: ``full``(All reactions were found), ``keyenzyme`` (Found key reactions indicate pathway presence (at least 66% of the pathway reactions are present)), ``NA`` (The pathway is not predicted to be present), ``threshold`` (Pathway is present because at least 80% of its reactions are present)
* ***For enzymes which have a high similar sequence to other enzymes with different function a higher identity cutoff is used for the blast search (this exceptions are defined in ``gapseq/dat/exceptions.tbl``)


## gapseq find-transport
### *-Transporter.tbl
Data about found transporter is listed in this tab-separated text file.

| id             | tc        | sub                   | exid        | rea                       | qseqid       | pident                          | evalue  | bitscore  | qcovs          | stitle        | sstart             | send             | 
| -              | -         | -                     | -           | -                         | -            | -                               | -       | -         | -              | -             | -                  | -                | 
| Transporter ID | TC number | Transported substance | Exchange ID | Associated model reaction | Query Seq-id | Percentage of identical matches | E value | Bit score | Query Coverage | Subject Title | Start of alignment | End of alignment |

## gapseq draft

### *-rxnWeights.RDS
Reaction weights table (temporary file needed for gapseq fill).

### *-rxnXgenes.RDS
Table with gene-X-reaction association (temporary file needed for gapseq fill).

### *-draft.RDS
Model draft file as R object.

### *-draft.xml
Draft model in SBML format.

## gapseq fill
### *.RDS
Final model saved as R object.
### *.xml
Final model in SBML format.
