addMetAttr <- function(mod, seed_x_mets) {
  require(data.table)
  require(stringr)
  
  #seed_x_mets <- fread(paste0(script.dir,"/../dat/seed_metabolites_edited.tsv"), header=T, stringsAsFactors = F, na.strings = c("null","","NA"))
  
  mettmp <- data.table(id = gsub("\\[.0\\]","", mod@met_id))
  
  mettmp <- merge(mettmp, seed_x_mets, sort = F, by = "id", all.x = T)
  
  #mettmp[, annotation := paste0("in_inchi;", inchikey)]
  
  colnames(mettmp)[which(colnames(mettmp)=="formula")] <- "chemicalFormula"
  
  # adding attributes
  mettmp$annotation <- "bqbiol_is;http://identifiers.org/SBO:0000247"
  
  # METANETX
  mettmp[!is.na(MNX_ID)  , annotation := paste0(annotation,";http://identifiers.org/metanetx.chemical/",MNX_ID)]
  
  # INCHI
  mettmp[!is.na(InChIKey), annotation := paste0(annotation,";http://identifiers.org/inchikey:",InChIKey)]
  mettmp[!is.na(InChI), annotation := paste0(annotation,";http://identifiers.org/inchi:",InChI)]
  
  # ModelSEED
  mettmp[!is.na(id) & !grepl("^cpd9", id), annotation := paste0(annotation,";http://identifiers.org/seed.compound/",id)]
  
  # HMDB
  tmp_ids <- str_split(mettmp[, hmdbID], pattern = ";")
  tmp_ids <- lapply(tmp_ids, FUN = function(x) {
    if(!is.na(x[1]))
      return(paste0(";http://identifiers.org/hmdb:",x, collapse = ""))
    else
      return(NA)
  })
  mettmp[, tmpids := tmp_ids]
  mettmp[!is.na(tmpids), annotation := paste0(annotation,tmpids)]
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
  mettmp[!is.na(tmpids), annotation := paste0(annotation,tmpids)]
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
  mettmp[!is.na(tmpids), annotation := paste0(annotation,tmpids)]
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
  mettmp[!is.na(tmpids), annotation := paste0(annotation,tmpids)]
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
  mettmp[!is.na(tmpids), annotation := paste0(annotation,tmpids)]
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
  mettmp[!is.na(tmpids), annotation := paste0(annotation,tmpids)]
  mettmp[, tmpids := NULL]; rm(tmp_ids)
  
  mod@met_attr <- as.data.frame(mettmp)
  
  return(mod)
}



