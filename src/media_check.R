library(data.table)
library(stringr)

media_check <- function(media.file, mod.orig, seed_x_mets){
  media   <- fread(media.file, stringsAsFactors = F, header = T, col.names = c("compounds", "name", "maxFlux"))
  bm.idx  <- grep("bio1", mod.orig@react_id)
  
  #
  # Check for very small numbers potentially causing issues with numeric accuracy
  #
  media.2small <- media[maxFlux<1e-6 & maxFlux > 0]
  if( nrow(media.2small) >0 ) warning(paste("Media components with very small numbers:", paste(media.2small$compunds, collapse = ",")))
  
  #
  # Check for presence of important biomass components (e.g. metals)
  #
  if( length(bm.idx) > 0 ){
    ess.met    <- c("cpd00009", "cpd00030", "cpd00058", "cpd00034", "cpd10515", "cpd10516","cpd00971","cpd00254", "cpd00205", "cpd00030", "cpd00149", "cpd00099", "cpd00063", "cpd01012", "cpd11574")
    bm.cpd     <- mod.orig@met_id[which(mod.orig@S[,bm.idx] != 0)]
    bm.ess.met <- intersect(ess.met, str_extract(bm.cpd, "cpd[0-9]+"))
    missing    <- ess.met[which(!bm.ess.met %in% media$compounds)]
    if( length(missing) > 0 ) warning(paste("Potentially missing compounds in medium:", paste(seed_x_mets[id %in% missing, paste(id, name)], collapse = ", ")))
  }
  
  #
  # Recommendations for cpd identifiers (degradation paths)
  #
  if( "cpd11601" %in% media$compounds ) warning("For pectin the identifier cpd27519 is recommended.")
  if( "cpd11657" %in% media$compounds ) warning("For starch the identifier cpd90003 is recommended.")
  
}