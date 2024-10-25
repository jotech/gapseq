add_missing_exchanges <- function(mod, ub = 1000) {
  ex.mets.ind  <- which(mod@met_comp==mod@mod_compart[2])
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
                    metComp = mod@mod_compart[2],
                    ub = ub,
                    lb = 0,
                    reactName = paste0(ex.mets.name[i], " Exchange"), 
                    metName = ex.mets.name[i])
    if("gs.origin" %in% colnames(mod@react_attr)) {
      mod@react_attr[which(mod@react_id == paste0("EX_",gsub("\\[.*\\]","",ex.mets.ids[i]),"_e0")),c("gs.origin","seed")] <- data.frame(gs.origin = 7,
                                                                                                                                        seed = paste0("EX_",gsub("\\[.*\\]","",ex.mets.ids[i]),"_e0"),
                                                                                                                                        stringsAsFactors = F)
    }
  }
  
  return(mod)
}

add_sinks <- function(mod, ub = 1000) {
  in.mets.ind  <- which(mod@met_comp==mod@mod_compart[1])
  in.mets.ids  <- mod@met_id[in.mets.ind]
  in.mets.name <- mod@met_name[in.mets.ind]
  
  for(i in 1:length(in.mets.ids)) {
    mod <- addReact(model = mod,
                    id = paste0("SINK_",gsub("\\[.*\\]","",in.mets.ids[i]),"_c0"), 
                    met = in.mets.ids[i],
                    Scoef = -1,
                    metComp = mod@mod_compart[1],
                    ub = ub,
                    lb = 0,
                    reactName = paste0(in.mets.name[i], " Sink"), 
                    metName = in.mets.name[i])
  }
  return(mod)
}

add_met_sink <- function(mod, cpd, obj = 0) {
  mname <- paste0(cpd,"_c0")
  if(paste0(cpd,"[c0]") %in% mod@met_id)
    mname <- NA_character_
  
  mod <- addReact(mod,
                  id = paste0("EX_",cpd,"_c0"),
                  met = paste0(cpd,"[c0]"),
                  Scoef = -1,
                  lb = 0,
                  ub = 1000,
                  reactName = paste0("EX ",cpd," c0"),
                  metName = mname,
                  metComp = mod@mod_compart[1],
                  obj = obj)
  return(mod)
}

add_missing_diffusion <- function(mod, ub = 1000){
  if( grepl("/src$|^src$", script.dir) ) {
    rxndb.path <- paste0(script.dir, "/../dat/diffusion_mets.tsv")
  }else{
    rxndb.path <- paste0(script.dir, "/dat/diffusion_mets.tsv")
  }
  #rxndb.path <- str_replace_all(rxndb.path, "/src","")
  diff.mets <- fread(rxndb.path, header=T, stringsAsFactors = F)
  mod <- add_reaction_from_db(mod, react = diff.mets$diffrxn, gs.origin = 8)
  mod <- add_missing_exchanges(mod)
  
  return(mod)
}

add_reaction_from_db <- function(mod, react, gs.origin = NA) {
  if( grepl("/src$|^src$", script.dir) ) {
    rxndb.path <- paste0(script.dir, "/../dat/seed_reactions_corrected.tsv")
  }else{
    rxndb.path <- paste0(script.dir, "/dat/seed_reactions_corrected.tsv")
  }

  #rxndb.path <- str_replace_all(rxndb.path, "/src","")
  mseed <- fread(rxndb.path, header=T, stringsAsFactors = F)
  mseed <- mseed[gapseq.status %in% c("approved","corrected")]
  mseed <- mseed[order(id)]
  mseed <- mseed[id %in% react]
  
  if(nrow(mseed)==0) {
    warning("None of the reactions ids are found in the database.")
    return(mod)
  }
  if(nrow(mseed) != length(react)){
    warning("Some of the specified reaction IDs are not found in database.") # TODO: Print wrong ids
  }
  
  for(i in (1:nrow(mseed))) {
    mets  <- unlist(str_split(string = mseed[i,compound_ids],pattern = ";"))
    rxn.info <- str_split(unlist(str_split(string = mseed[i,stoichiometry],pattern = ";")), pattern = ":", simplify = T)
    
    met.comp  <- rxn.info[,3]
    met.comp.n <- ifelse(met.comp==0,"c0","e0")
    met.comp.n <- ifelse(met.comp>=2,"p0",met.comp.n)
    
    met.ids   <- paste0(rxn.info[,2],"[",met.comp.n,"]")
    met.scoef <- as.numeric(rxn.info[,1])
    met.name  <- str_replace_all(rxn.info[,5],"\\\"","")
    
    met.name  <- paste0(met.name,"-",met.comp.n)
    
    ind <- which(met.name=="" | is.na(met.name))
    met.name[ind] <- met.ids[ind]
    
    is.rev <- ifelse(mseed[i,reversibility] %in% c("<","="),T,F)
    only.backwards <- ifelse(mseed[i,reversibility]=="<",T,F)
    
    # ind.new.mets <- which(met.ids %in% mod@met_id)
    # ind.old.mets <- which(mod@met_id %in% met.ids[ind.new.mets])
    # 
    # met.name[ind.new.mets] <- mod@met_name[ind.old.mets]
    
    ind.notnew.mets <- which(met.ids %in% mod@met_id)
    if(length(ind.notnew.mets) > 0)
      met.name[ind.notnew.mets] <- NA_character_
    
    new.react <- !any(grepl(paste0(mseed[i,id],"_c0"),mod@react_id))
    
    mod <- addReact(model = mod, 
                    id = paste0(mseed[i,id],"_c0"), 
                    met = met.ids,
                    Scoef = met.scoef,
                    metComp = mod@mod_compart[as.integer(met.comp)+1],
                    ub = ifelse(only.backwards, 0, 1000),
                    lb = ifelse(is.rev, -1000, 0),
                    reactName = mseed[i, name], 
                    metName = met.name)
    if(new.react) {
      mod@react_attr[which(mod@react_id == paste0(mseed[i,id],"_c0")),c("gs.origin","seed")] <- data.frame(gs.origin = gs.origin,
                                                                                                           seed = mseed[i,id],
                                                                                                           stringsAsFactors = F)
    }
  }
  return(mod)
}

# do I still use this?
add_exchanges <- function(mod, cpd, ub = 1000, metname=NA) {
  metname <- gsub("-e0$","",metname)
  metname <- ifelse(!is.na(metname),paste0(metname,"-e0"),metname)
  
  for(m in cpd){
    ex.id  <- paste0("EX_",gsub("\\[.*\\]","",m),"_e0")
    if( ex.id %in% mod@react_id )
      next
    ex.met.id    <- paste0(gsub("\\[.*\\]","",m), "[e0]")
    ex.mets.ind  <- grep(m, mod@met_id,fixed=T)[1]
    if( all(is.na(metname)) )
      ex.mets.name <- mod@met_name[ex.mets.ind]
    else
      ex.mets.name <- metname[match(m, cpd)]
    mod <- addReact(model = mod,
                    id = ex.id, 
                    met = ex.met.id,
                    Scoef = -1,
                    metComp = mod@mod_compart[2],
                    ub = ub,
                    lb = 0,
                    reactName = paste0(ex.mets.name, " Exchange"), 
                    metName = ex.mets.name)
    mod@react_attr[which(mod@react_id == ex.id), c("gs.origin","seed")] <- data.frame(gs.origin = 7,
                                                                                      seed = ex.id,
                                                                                      stringsAsFactors = F)
  }
  
  return(mod)
}

# find orphan exchange reactions that are not linked to other reactions
rm_unused_exchanges <- function(mod) {
  mod <- copy(mod)
  present.exchanges <- str_match(mod@react_id[grep("^EX_.*_e0",mod@react_id)],"_(.*?)_")
  ex.using.rxn <- sapply(present.exchanges[,2], function(exmet){length(mod@react_id[which(mod@S[match(paste0(exmet,"[e0]"),mod@met_id),]!=0)])})
  using.count <- unname(which(ex.using.rxn < 2))
  if(identical(using.count, integer(0))) return(mod)
  ex.unused <- paste0("EX_", present.exchanges[,2][using.count],"_e0") # external metabolites only used by one reactions
  if(!all(ex.unused %in% mod@react_id)) warning("Exchange reactions not found")
  mod <- rmReact(mod, react=ex.unused)
  return(mod)
}
