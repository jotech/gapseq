addGeneAttr <- function(mod, dtg) {
  dtg <- dtg[order(stitle, -bitscore)][!duplicated(stitle)]
  dtg <- dtg[,.(qseqid, stitle)]
  
  dtg$cvt <- NA_character_
  
  # UniRef
  dtg[grepl("^UniRef(100|90|50)", qseqid), cvt := paste0("uniref:",qseqid)]
  dtg[grepl("^UniRef(100|90|50)", qseqid), cvt := gsub("\\|.*$","",cvt)]
  
  # Uniprot from trembl/swissprot
  dtg[grepl("^tr\\||^sp\\|", qseqid), cvt := paste0("uniprot:",unlist(strsplit(qseqid,"\\|"))[2]), by = stitle]
  
  # Uniprot from TC-DB
  dtg[grepl("\\|TC-DB\\|",qseqid), cvt := paste0("uniprot:",unlist(strsplit(qseqid,"\\|"))[3]), by = stitle]
  
  # chekc if uniref IDs are correct
  dtg[grepl("^uniref",cvt)  & !grepl("UniRef(100|90|50)_([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}|UPI[A-F0-9]{10})$",cvt), cvt := NA_character_]
  
  # check if uniprot keys are correct
  dtg[grepl("^uniprot",cvt) & !grepl("^uniport\\:([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\\.\\d+)?$",cvt), cvt := NA_character_]
  
  inds <- match(mod@allGenes, dtg$stitle)
  
  mod@genes_attr$CVTerms <- ifelse(!is.na(dtg$cvt[inds]) & !grepl("bqbiol_isHomologTo",dtg$cvt[inds]),
                                   paste0("bqbiol_isHomologTo;http://identifiers.org/",dtg$cvt[inds]),
                                   NA_character_)
  
  return(mod)
}

