# What kind of organisms are supported by gapseq?
By default gapseq searches for pathways using bacterial sequences.
This can be changed to other taxonomic ranks, for example to animals `gapseq find -t Metazoa`.
The taxonomic rank should be compatible with [uniprot](https://www.uniprot.org/help/taxonomy) because the respective sequences are download from uniprot.
However, the reconstruction and gapfilling of metabolic networks is currently limited to bacteria and archaea.

# What are the output files of gapseq?
Predicted pathways, reactions, and transporters are stored in tabulator-separated text files (.tbl) and the final metabolic network is saved as [systems biology markup language (SBML)](http://sbml.org) (.xml), which is an xml-based markup language for biological model exchange.
In addition, final output and temporary files are stored in the [R Data Format (RDS)](http://www.sthda.com/english/wiki/saving-data-into-r-data-format-rds-and-rdata).

# What are the RDS output files for?
Temporary model files are stored in [R Data Format (RDS)](http://www.sthda.com/english/wiki/saving-data-into-r-data-format-rds-and-rdata) (files containing the R-object 'ModelOrg' from the R-package 'cobrar') and serve as input for other substeps such as gapseq fill and adapt. Since users can do further analysis with R, e.g. using [cobrar](https://github.com/Waschina/cobrar), the final output files are also in RDS format.

# How is pathway completeness calculated in gapeq?

The `gapseq find` output file "[...]-Pathways.tbl" contains a column named 'Completeness', which states the predicted pathway completeness in percent. The completeness is calculated using the formula:

$$
C = \frac{N_{found}}{N_{pwy} - N_{vague} - N_{spont}}
$$

Where, :math:`N_{found}` is the number of reactions found, :math:`N_{pwy}` is the total number of reactions in the pathway, :math:`N_{vague}` is the number of reactions without reference sequences, and $N_{spont}$ is the number of spontaneous reactions in the pathway.

In cases, where the number of vague reactions (:math:`N_{vague}`) accounts for 30% or more of all pathway reactions, the formula is modified to:

$$
C = \frac{N_{found}}{N_{pwy} - N_{spont}}
$$

This prevents over-prediction of pathways where the majority of reactions cannot be predicted due to missing reference sequences.

# Which codon table does gapseq use when translating a nucleotide genome?

gapseq automatically selects the appropriate codon translation table by running pyrodigal with three options:

  - Translation Table 4: ["Mycoplasma/Spiroplasma (Mollicutes)"](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes#SG4)
  - Translation Table 11: ["Bacterial, Archaeal, and Plant Plastid Code"](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes#SG11) (default for most prokaryotic tools)
  - Translation Table 25: ["Candidate Division SR1 and Gracilibacteria"](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes#SG25)

The choice between Table 11 and Tables 4/25 depends on genome coverage. If using Table 4 or 25 gives at least 5% higher coverage than Table 11, then 4 or 25 is used. Choosing between Table 4 and 25 is more nuanced since both yield the same coverage. The key difference is how the codon UGA is interpreted:

  - In Table 11, UGA is a stop codon.
  - In Table 4, UGA codes for Tryptophan.
  - In Table 25, UGA codes for Glycine.

Since the Tryptophan content in proteins is typically around 1%, the table that produces a Tryptophan usage closest to this value is selected.

Admittedly, this approach relies on heuristic thresholds, but it works well in practice. If users already know the correct codon table for their genome, they can provide a protein FASTA file directly to avoid translation by gapseq.
