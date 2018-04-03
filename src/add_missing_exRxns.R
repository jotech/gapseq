add_missing_exchanges <- function(mod, ub = 1000) {
  mod <- copy(mod)
  ex.mets.ind  <- which(mod@met_comp==2)
  ex.mets.ids  <- mod@met_id[ex.mets.ind]
  ex.mets.name <- mod@met_name[ex.mets.ind]
  
  present.exchanges <- str_match(mod@react_id[grep("^EX_.*_e0",mod@react_id)],"_(.*?)_")
  if(ncol(present.exchanges)==2){
    ind.new <- !(gsub("\\[.*\\]","",ex.mets.ids) %in% present.exchanges[,2])
    ex.mets.ind  <- ex.mets.ind[ind.new]
    ex.mets.ids  <- ex.mets.ids[ind.new]
    ex.mets.name <- ex.mets.name[ind.new]
  }

  if(length(ex.mets.ids)==0) return(mod)
  
  for(i in 1:length(ex.mets.ids)) {
    #cat("Adding new exchange reaction for: ",ex.mets.ids[i],"\n")
    mod <- addReact(model = mod,
                    id = paste0("EX_",gsub("\\[.*\\]","",ex.mets.ids[i]),"_e0"), 
                    met = ex.mets.ids[i],
                    Scoef = -1,
                    reversible = T, 
                    metComp = 2,
                    ub = ub,
                    lb = 0,
                    reactName = paste0(ex.mets.name[i], " Exchange"), 
                    metName = ex.mets.name[i])
  }
  
  return(mod)
}

add_sinks <- function(mod, ub = 1000) {
  in.mets.ind  <- which(mod@met_comp==1)
  in.mets.ids  <- mod@met_id[in.mets.ind]
  in.mets.name <- mod@met_name[in.mets.ind]
  
  for(i in 1:length(in.mets.ids)) {
    mod <- addReact(model = mod,
                    id = paste0("SINK_",gsub("\\[.*\\]","",in.mets.ids[i]),"_c0"), 
                    met = in.mets.ids[i],
                    Scoef = -1,
                    reversible = F, 
                    metComp = 1,
                    ub = ub,
                    lb = 0,
                    reactName = paste0(in.mets.name[i], " Sink"), 
                    metName = in.mets.name[i])
  }
  return(mod)
}

add_met_sink <- function(mod, cpd, obj = 0) {
  mod <- addReact(mod,
                  id = paste0("EX_",cpd,"_c0"),
                  met = paste0(cpd,"[c0]"),
                  Scoef = -1,
                  reversible = F,
                  lb = 0,
                  ub = 1000,
                  reactName = paste0("EX ",cpd," c0"),
                  metName = paste0(cpd,"_c0"),
                  metComp = 1,
                  obj = obj)
  return(mod)
}