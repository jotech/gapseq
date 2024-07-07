constrain.model <- function(mod, media.file=NA, media=NA, scaling.fac = 1, ub = 1000) {
  if( all(is.na(media)) & all(is.na(media.file)))
    stop("Needed media file or media data frame")
  if( !all(is.na(media)) & !all(is.na(media.file)))
    warning("Media file and media data frame provided. Data frame will be used.")
  if( all(is.na(media)) )
    media <- fread(media.file, stringsAsFactors = F, header = T, col.names = c("compounds", "name", "maxFlux"))
  if( is.na(media.file) ) media <- copy(media) # copy data.table before changing it
  
  ex.rxns <- grep("^EX_",mod@react_id, fixed = F)
  mod@lowbnd[ex.rxns] <- 0
  media[,compounds := paste0("EX_",compounds,"_e0")]
  media$mod.rxn.id <- match(media$compounds,mod@react_id)
  
  # Metabolited alread present
  media1 <- copy(media[!is.na(mod.rxn.id)])
  mod@lowbnd[media1$mod.rxn.id] <- -media1$maxFlux * scaling.fac

  # New metabolites (new exchange reactions need to be added)
  media2 <- copy(media[is.na(mod.rxn.id)])
  media2[,compounds2 := gsub("EX_","",compounds)]
  media2[,compounds2 := gsub("_.0","",compounds2)]
  gs.origin.tag <- !is.null(mod@react_attr$gs.origin) # If this tag exists in reaction attribute table, use it!
  if(nrow(media2)>0) {
    for(i in 1:nrow(media2)) {
      mod <- addReact(model = mod,
                      id = media2[i,compounds],
                      met = paste0(media2[i,compounds2],"[e0]"),
                      Scoef = -1,
                      metComp = mod@mod_compart[2],
                      ub = ub,
                      lb = -media2[i, maxFlux] * scaling.fac,
                      reactName = paste0(media2[i,name], " Exchange"), 
                      metName = media2[i,name])
      if(gs.origin.tag)
        mod@react_attr[which(mod@react_id == media2[i,compounds]),c("gs.origin","seed")] <- c(7,gsub("_.0","",media2[i,compounds]))
      
    }
  }
  return(mod)
}
