addReactAttr <- function(mod) {
  require(data.table)
  require(stringr)

  seed_x_metCyc <- fread(paste0(script.dir,"/../dat/mnxref_seed-other.tsv"), header = T)
    
  reacttmp <- data.table(mod@react_attr)
  reacttmp[is.na(CVTerms), CVTerms := ""]
  reacttmp$seed <- mod@react_id
  reacttmp[grepl("^rxn",seed), seed := gsub("_.0$","",seed)]
  #print(reacttmp)
  
  
  # # adding general attributes
  # reacttmp$annotation <- "bqbiol_is;http://identifiers.org/SBO:0000167" # "biochemical or transport reaction"
  
  # EC-code
  ec_id_mult <- function(ec) {
    ecs <- unique(unlist(strsplit(ec, "/", fixed = T)))
    ecs <- ecs[ecs!=""]
    
    ecs_out <- paste0("http://identifiers.org/ec-code:",ecs, collapse = ";")
    
    return(ecs_out)
  }
  reacttmp[!is.na(ec) & ec!="", CVTerms := paste0(CVTerms,";",ec_id_mult(ec)), by = seed]
  
  # SEED
  reacttmp[!is.na(seed) & seed!="" & grepl("^rxn", seed) & !grepl("^rxn9", seed), 
           CVTerms := paste0(CVTerms,";http://identifiers.org/seed.reaction:",seed)]
  
  # # SBO Term for exchange reactions 
  # reacttmp[grepl("^EX_", seed), annotation := paste0(annotation,";http://identifiers.org/SBO:0000627")]
  
  # # demand reactions
  # reacttmp[grepl("^DM_",mod@react_id), annotation := paste0(annotation,";http://identifiers.org/SBO:0000628")]
  # reacttmp[grepl("^rxn13782",mod@react_id), annotation := paste0(annotation,";http://identifiers.org/SBO:0000628")]
  # reacttmp[grepl("^rxn13783",mod@react_id), annotation := paste0(annotation,";http://identifiers.org/SBO:0000628")]
  # reacttmp[grepl("^rxn13784",mod@react_id), annotation := paste0(annotation,";http://identifiers.org/SBO:0000628")]
  
  # # annotate transport reactions
  # prod <- apply(mod@S, 2, FUN = function(x) which(x > 0))
  # cons <- apply(mod@S, 2, FUN = function(x) which(x < 0))
  # prod <- lapply(prod, FUN = function(x) sort(gsub("\\[.0\\]","",mod@met_id[x])))
  # cons <- lapply(cons, FUN = function(x) sort(gsub("\\[.0\\]","",mod@met_id[x])))
  # prod <- lapply(prod, FUN = function(x) paste(x, collapse = "$"))
  # cons <- lapply(cons, FUN = function(x) paste(x, collapse = "$"))
  # is.transp <- which(unlist(lapply(1:react_num(mod), FUN = function(x) prod[[x]] == cons[[x]])))
  # reacttmp[is.transp, annotation := paste0(annotation,";http://identifiers.org/SBO:0000185")]
  # # TODO: Find ATPase transporters
  
  # # anything else: biochemical reaction (SBO:0000176)
  # reacttmp[!grepl("SBO:0000185",annotation) & !grepl("SBO:0000627",annotation), 
  #          annotation := paste0(annotation,";http://identifiers.org/SBO:0000176")]
  
  
  #~~~~~~~~~~~~~~~#
  # X-referencing #
  #~~~~~~~~~~~~~~~#
  mx_reactxref <- fread(paste0(paste0(script.dir,"/../dat/mnxref_reac_xref.tsv")), skip = 386)
  colnames(mx_reactxref) <- c("XREF", "MNX_ID","sr")
  mx_reactxref[, sr := NULL]
  
  # METANETX
  mx_seed <- copy(seed_x_metCyc)
  mx_seed[, other := NULL]
  mx_seed <- mx_seed[!duplicated(paste(MNX_ID,seed, sep = "$"))]
  mx_seed <- mx_seed[!duplicated(mx_seed)]
  mx_seed[, seedID := seed]
  
  # KEGG
  mx_kegg <- copy(mx_reactxref[grepl("^kegg:",XREF)])
  mx_kegg[, XREF := gsub("^kegg:","", XREF)]
  mx_kegg <- mx_kegg[grepl("^R[0-9]{5}", XREF)]
  mx_kegg <- mx_kegg[, .(keggID = paste(XREF, collapse = ";")), by = MNX_ID]
  mx_seed <- merge(mx_seed, mx_kegg, by = "MNX_ID", all.x = T, sort = F)
  
  # bigg
  mx_bigg <- copy(mx_reactxref[grepl("^bigg:",XREF)])
  mx_bigg[, XREF := gsub("^bigg:","", XREF)]
  mx_bigg <- mx_bigg[!grepl("^R_", XREF)]
  mx_bigg <- mx_bigg[!grepl("^EX_", XREF)] 
  mx_bigg <- mx_bigg[, .(biggID = paste(XREF, collapse = ";")), by = MNX_ID]
  mx_seed <- merge(mx_seed, mx_bigg, by = "MNX_ID", all.x = T, sort = F)
  
  # biocyc
  mx_biocyc <- copy(mx_reactxref[grepl("^metacyc:",XREF)])
  mx_biocyc[, XREF := gsub("^metacyc:","", XREF)]
  mx_biocyc <- mx_biocyc[, .(biocycID = paste(XREF, collapse = ";")), by = MNX_ID]
  mx_seed   <- merge(mx_seed, mx_biocyc, by = "MNX_ID", all.x = T, sort = F)
  
  # delete old mappings if they exist
  if("seedID" %in% colnames(reacttmp)) {
    reacttmp[, biocycID := NULL]; reacttmp[, keggID := NULL]; reacttmp[, biggID := NULL]
    reacttmp[, seedID := NULL]; reacttmp[, MNX_ID := NULL]
  }
  
  # add X-refs to annotations
  reacttmp <- merge(reacttmp, mx_seed, by = "seed", sort = F, all.x = T)
  # ModelSEED
  reacttmp[!is.na(seedID) & !grepl("^rxn9", seedID), CVTerms := paste0(CVTerms,";http://identifiers.org/seed.reaction:",seedID)]
  # kegg
  tmp_ids <- str_split(reacttmp[, keggID], pattern = ";")
  tmp_ids <- lapply(tmp_ids, FUN = function(x) {
    if(!is.na(x[1]))
      return(paste0(";http://identifiers.org/kegg.reaction:",x, collapse = ""))
    else
      return(NA)
  })
  reacttmp[, tmpids := tmp_ids]
  reacttmp[!is.na(tmpids), CVTerms := paste0(CVTerms,tmpids)]
  reacttmp[, tmpids := NULL]; rm(tmp_ids)
  # bigg
  tmp_ids <- str_split(reacttmp[, biggID], pattern = ";")
  tmp_ids <- lapply(tmp_ids, FUN = function(x) {
    if(!is.na(x[1]))
      return(paste0(";http://identifiers.org/bigg.reaction:",x, collapse = ""))
    else
      return(NA)
  })
  reacttmp[, tmpids := tmp_ids]
  reacttmp[!is.na(tmpids), CVTerms := paste0(CVTerms,tmpids)]
  reacttmp[, tmpids := NULL]; rm(tmp_ids)
  # biocyc
  tmp_ids <- str_split(reacttmp[, biocycID], pattern = ";")
  tmp_ids <- lapply(tmp_ids, FUN = function(x) {
    if(!is.na(x[1]))
      return(paste0(";http://identifiers.org/biocyc:META:",x, collapse = ""))
    else
      return(NA)
  })
  reacttmp[, tmpids := tmp_ids]
  reacttmp[!is.na(tmpids), CVTerms := paste0(CVTerms,tmpids)]
  reacttmp[, tmpids := NULL]; rm(tmp_ids)
  # metanetx
  tmp_ids <- str_split(reacttmp[, MNX_ID], pattern = ";")
  tmp_ids <- lapply(tmp_ids, FUN = function(x) {
    if(!is.na(x[1]))
      return(paste0(";http://identifiers.org/metanetx.reaction:",x, collapse = ""))
    else
      return(NA)
  })
  reacttmp[, tmpids := tmp_ids]
  reacttmp[!is.na(tmpids), CVTerms := paste0(CVTerms,tmpids)]
  reacttmp[, tmpids := NULL]; rm(tmp_ids)
  
  
  # Export
  reacttmp[!is.na(CVTerms) & CVTerms != "" & !grepl("^bqbiol_is", CVTerms), CVTerms := paste0("bqbiol_is",CVTerms)]
  mod@react_attr <- as.data.frame(reacttmp)
  
  return(mod)
}
