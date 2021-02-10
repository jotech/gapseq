# Unit tests

This folder contains script to check basic and advanced functions of gapseq.

## model_flux
The script ``model_flux.sh`` checks if model predictions meet certain expected phenotypes as defined in ``model_flux.tbl``.
For every organism, the genome is downloaded (if necessary) to the folder ``genomes/`` and the genome serves as input for gapseq to create a metabolic model which is saved, together with all intermediate files to ``out/``.
Finally, the final model is used to make flux predictions that are compared to expected phenotypic data as defined in ``model_flux.tbl``.
All results are appended to the file model_flux.log including the gapseq version number in order to track the performance of different program releases.

