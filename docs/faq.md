# What kind of organisms are supported by gapseq?
By default gapseq searches for pathways using bacterial sequences.
This can be changed to other taxonomic ranks, for example to animals `gapseq find -t Metazoa`.
The taxonomic rank should be compatible with [uniprot](https://www.uniprot.org/help/taxonomy) because the respective sequences are download from uniprot.
However, the reconstruction and gapfilling of metabolic networks is currently limited to bacteria and archaea.

# What are the output files of gapseq?
Predicted pathways, reactions, and transporters are stored in tabulator-separated text files (.tbl) and the final metabolic network is saved as [systems biology markup language (SBML)](http://sbml.org) (.xml), which is an xml-based markup language for biological model exchange.
In addition, final output and temporary files are stored in the [R Data Format (RDS)](http://www.sthda.com/english/wiki/saving-data-into-r-data-format-rds-and-rdata).

# What are all these RDS output files for?
Temporary files are stored in [R Data Format (RDS)](http://www.sthda.com/english/wiki/saving-data-into-r-data-format-rds-and-rdata) and serve as input for other substeps such as gapseq draft, fill, or adapt.
Since we propagate the further analysis with R using for example [cobrar](https://github.com/Waschina/cobrar), the final output files are also kept in RDS format.

# How is pathway completeness calculated in gapeq?

The `gapseq find` output file "[...]-Pathways.tbl" contains a column named 'Completeness', which states the predicted pathway completeness in percent. The completeness is calculated using the formula:

$$
C = \frac{N_{found}}{N_{pwy} - N_{vague} - N_{spont}}
$$
Where, $N_{found}$ is the number of reactions found, $N_{pwy}$ is the total number of reactions in the pathway, $N_{vague}$ is the number of reactions without reference sequences, and $N_{spont}$ is the number of spontaneous reactions in the pathway.

In cases, where the number of vague reactions ($N_{vague}$) accounts for 30% or more of all pathway reactions, the formula is modified to:
$$
C = \frac{N_{found}}{N_{pwy} - N_{spont}}
$$
This prevents over-prediction of pathways where the majority of reactions cannot be predicted due to missing reference sequences.
