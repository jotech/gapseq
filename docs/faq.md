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

