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


add_missing_diffusion <- function(mod, ub = 1000){
  diff.mets.ids <- c("cpd00011", #"co2
                     "cpd00007", #o2
                     "cpd11640", #"h2
                     "cpd00001", #h2o
                     "cpd00363", # etoh
                     "cpd00055", # Formaldehyde
                     "cpd01861", # (R)-1,2-Propanediol
                     "cpd00453", # 1,2-Propanediol
                     "cpd00150", # HCN (Cyanide)
                     "cpd00448", # D-Glyceraldehyde
                     "cpd00025", # H2O2
                     "cpd00359", # Indole (10.1128/JB.01477-10 <https://dx.doi.org/10.1128%2FJB.01477-10>)
                     "cpd00811", # cpd00811 / tmao (klein und ungeladen)
                     "cpd12733", # SO2 (klein, ungeladen)
                     "cpd00371", # Propanal (klein, ungeladen)
                     "cpd00418", # NO
                     "cpd00659", # N2O
                     "cpd11574", # Molybdate  (H2MoO4)
                     "cpd00448", # glyald (D-Glyceraldehyde)
                     "cpd00450", # dms (Methyl sulfide)
                     "cpd08021", # DMSO
                     "cpd00071", # acald (Acetaldehyde)
                     "cpd00281", # 4abut (GABA)
                     "cpd03828" # 23dappa (2,3-Diaminopropionate))
                    )
  diff.mets.ids.comp <- paste0(diff.mets.ids,"[c0]")
  diff.mets.names <- ifelse(diff.mets.ids.comp %in% mod@met_id,
                            mod@met_name[match(diff.mets.ids.comp, mod@met_id)],
                            diff.mets.ids.comp)
  
  for(i in 1:length(diff.mets.ids)) {
    diff.ex.id <- paste0("EX_",gsub("\\[.*\\]","",diff.mets.ids[i]),"_e0")
    if( diff.ex.id %in% mod@react_id)
      mod <- rmReact(mod, react=diff.ex.id)
      
    mod <- addReact(model = mod,
                    id = diff.ex.id, 
                    met = diff.mets.ids.comp[i],
                    Scoef = -1,
                    reversible = T, 
                    metComp = 2,
                    ub = ub,
                    lb = 0,
                    reactName = paste0(diff.mets.names[i], " Exchange + Diffusion"), 
                    metName = diff.mets.names[i])
  }
  return(mod)
}