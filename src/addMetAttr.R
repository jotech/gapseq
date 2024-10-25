addMetAttr <- function(mod, seed_x_mets) {
  require(data.table)
  require(stringr)
  
  #seed_x_mets <- fread(paste0(script.dir,"/../dat/seed_metabolites_edited.tsv"), header=T, stringsAsFactors = F, na.strings = c("null","","NA"))
  
  # mettmp <- data.table(id = gsub("\\[.0\\]","", mod@met_id))
  mettmp <- data.table(mod@met_attr)[,.(CVTerms, SBOTerm)]
  mettmp$id <- gsub("\\[.0\\]","", mod@met_id)
  
  mettmp <- merge(mettmp, seed_x_mets, sort = F, by = "id", all.x = T)
  setnames(mettmp,"formula","chemicalFormula")
  
  # adding attributes
  mettmp[is.na(CVTerms), CVTerms := ""]
  
  # METANETX
  mettmp[!is.na(MNX_ID)  , CVTerms := paste0(CVTerms,";http://identifiers.org/metanetx.chemical/",MNX_ID)]
  
  # INCHI
  mettmp[!is.na(InChIKey), CVTerms := paste0(CVTerms,";http://identifiers.org/inchikey:",InChIKey)]
  
  # ModelSEED
  mettmp[!is.na(id) & !grepl("^cpd9", id), CVTerms := paste0(CVTerms,";http://identifiers.org/seed.compound/",id)]
  
  # HMDB
  tmp_ids <- str_split(mettmp[, hmdbID], pattern = ";")
  tmp_ids <- lapply(tmp_ids, FUN = function(x) {
    if(!is.na(x[1]))
      return(paste0(";http://identifiers.org/hmdb:",x, collapse = ""))
    else
      return(NA)
  })
  mettmp[, tmpids := tmp_ids]
  mettmp[!is.na(tmpids), CVTerms := paste0(CVTerms,tmpids)]
  mettmp[, tmpids := NULL]; rm(tmp_ids)
  
  # reactome
  tmp_ids <- str_split(mettmp[, reactomeID], pattern = ";")
  tmp_ids <- lapply(tmp_ids, FUN = function(x) {
    if(!is.na(x[1]))
      return(paste0(";http://identifiers.org/reactome:R-ALL-",x, collapse = ""))
    else
      return(NA)
  })
  mettmp[, tmpids := tmp_ids]
  mettmp[!is.na(tmpids), CVTerms := paste0(CVTerms,tmpids)]
  mettmp[, tmpids := NULL]; rm(tmp_ids)
  
  # kegg
  tmp_ids <- str_split(mettmp[, keggID], pattern = ";")
  tmp_ids <- lapply(tmp_ids, FUN = function(x) {
    if(!is.na(x[1]))
      return(paste0(";http://identifiers.org/kegg.compound:",x, collapse = ""))
    else
      return(NA)
  })
  mettmp[, tmpids := tmp_ids]
  mettmp[!is.na(tmpids), CVTerms := paste0(CVTerms,tmpids)]
  mettmp[, tmpids := NULL]; rm(tmp_ids)
  
  # chebi
  tmp_ids <- str_split(mettmp[, chebiID], pattern = ";")
  tmp_ids <- lapply(tmp_ids, FUN = function(x) {
    if(!is.na(x[1]))
      return(paste0(";http://identifiers.org/CHEBI:",x, collapse = ""))
    else
      return(NA)
  })
  mettmp[, tmpids := tmp_ids]
  mettmp[!is.na(tmpids), CVTerms := paste0(CVTerms,tmpids)]
  mettmp[, tmpids := NULL]; rm(tmp_ids)
  
  # bigg
  tmp_ids <- str_split(mettmp[, biggID], pattern = ";")
  tmp_ids <- lapply(tmp_ids, FUN = function(x) {
    if(!is.na(x[1]))
      return(paste0(";http://identifiers.org/bigg.metabolite:",x, collapse = ""))
    else
      return(NA)
  })
  mettmp[, tmpids := tmp_ids]
  mettmp[!is.na(tmpids), CVTerms := paste0(CVTerms,tmpids)]
  mettmp[, tmpids := NULL]; rm(tmp_ids)
  
  # biocyc
  tmp_ids <- str_split(mettmp[, biocycID], pattern = ";")
  tmp_ids <- lapply(tmp_ids, FUN = function(x) {
    if(!is.na(x[1]))
      return(paste0(";http://identifiers.org/biocyc:META:",x, collapse = ""))
    else
      return(NA)
  })
  mettmp[, tmpids := tmp_ids]
  mettmp[!is.na(tmpids), CVTerms := paste0(CVTerms,tmpids)]
  mettmp[, tmpids := NULL]; rm(tmp_ids)
  
  mettmp[!is.na(CVTerms) & CVTerms != "" & !grepl("^bqbiol_is", CVTerms), CVTerms := paste0("bqbiol_is",CVTerms)]
  mod@met_attr <- as.data.frame(mettmp)
  
  return(mod)
}
