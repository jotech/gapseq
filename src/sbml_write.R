write_gapseq_sbml <- function(mod, out.id, verbose=FALSE, overwrite=TRUE) {
  
  # handling of missing values/attributes
  mod@met_attr$chemicalFormula <- ifelse(mod@met_attr$chemicalFormula=="null",
                                         "",
                                         mod@met_attr$chemicalFormula)
  
  # handling of subsystem names (and remove empty subsystem)
  colnames(mod@subSys) <- gsub("^\\||\\|$","",colnames(mod@subSys))
  
  # Model notes to include gapseq version
  notes_str <- paste0("<notes>\n  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>",mod@mod_desc,"</p>\n  </html>\n</notes>")
  # mod@mod_attr <- data.frame(notes = notes_str)
  mod@mod_notes <- notes_str
  
  # writing sbml
  sbml.f <- paste0(out.id, ".xml")
  if(overwrite==FALSE & file.exists(sbml.f)){
    out.tmp <- tools::file_path_sans_ext(sbml.f)
    out.tmp <- make.unique(c(out.tmp, out.tmp), sep="")[2] # add increasing number to existing output file
    sbml.f <- paste0(out.tmp, ".xml")
  }
  if(verbose) print(paste0("Writing SBML file ", sbml.f))
  sbml.o <- cobrar::writeSBMLmod(mod, file_path = sbml.f)
  if(sbml.o==F)
    warning("Writing SBML-file failed.")
}
