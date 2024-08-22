write_gapseq_sbml <- function(mod, out.id) {
  
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
  sbml.o <- cobrar::writeSBMLmod(mod, file_path = paste0(out.id, ".xml"))
  if(sbml.o==F)
    warning("Writing SBML-file failed.")
}
