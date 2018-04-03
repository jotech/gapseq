constrain.model <- function(mod, media.file, scaling.fac = 1, ub = 1000) {
  ex.rxns <- grep("^EX_",mod@react_id, fixed = F)
  mod@lowbnd[ex.rxns] <- 0
  media <- fread(media.file, stringsAsFactors = F, header = T)
  media[,compounds := paste0("EX_",compounds,"_e0")]
  media$mod.rxn.id <- match(media$compounds,mod@react_id)
  
  # Metabolited alread present
  media1 <- copy(media[!is.na(mod.rxn.id)])
  mod@lowbnd[media1$mod.rxn.id] <- -media1$maxFlux * scaling.fac

  # New metabolites (new exchange reactions need to be added)
  media2 <- copy(media[is.na(mod.rxn.id)])
  media2[,compounds2 := gsub("EX_","",compounds)]
  media2[,compounds2 := gsub("_.0","",compounds2)]
  if(nrow(media2)>0) {
    for(i in 1:nrow(media2)) {
      mod <- addReact(model = mod,
                      id = media2[i,compounds], 
                      met = paste0(media2[i,compounds2],"[e0]"),
                      Scoef = -1,
                      reversible = T, 
                      metComp = 2,
                      ub = ub,
                      lb = -media2[i, maxFlux] * scaling.fac,
                      reactName = paste0(media2[i,name], " Exchange"), 
                      metName = media2[i,name])
    }
  }
  return(mod)
}