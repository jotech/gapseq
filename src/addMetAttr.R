addMetAttr <- function(mod, seed_x_mets) {
  require(data.table)
  require(stringr)
  
  #seed_x_mets <- fread(paste0(script.dir,"/../dat/seed_metabolites_edited.tsv"), header=T, stringsAsFactors = F, na.strings = c("null","","NA"))
  
  mettmp <- data.table(id = gsub("\\[.0\\]","", mod@met_id))
  
  mettmp <- merge(mettmp, seed_x_mets, sort = F, by = "id", all.x = T)
  
  mod@met_attr <- as.data.frame(mettmp)
  
  return(mod)
}



