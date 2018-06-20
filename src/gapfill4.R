gapfill4 <- function(mod.orig, mod.full, core.rxn, min.gr = 0.1, dummy.bnd = 1e-3, diet.scale = 1,
                     core.weight = 1, dummy.weight = 10000, script.dir, core.only = FALSE, 
                     mtf.scale = 3, verbose=verbose) {
  source(paste0(script.dir, "/src/sysBiolAlg_mtfClass2.R"))
  # backup model
  mod.orig.bak <- mod.orig
  
  # Load dummy model
  mod <- mod.full
  mod <- sync_full_mod(mod.orig, mod)
  
  # get currently present reactions
  pres.rxns <- mod.orig@react_id
  pres.rxns <- gsub("_.*","",pres.rxns)
  
  # Get all reactions (non-duplicates) that are not yet part of the model
  mseed <- fread(paste0(script.dir, "/dat/seed_reactions_corrected.tsv"), header=T, stringsAsFactors = F)
  mseed <- mseed[gapseq.status %in% c("approved","corrected")]
  mseed <- mseed[!(id %in% pres.rxns)]
  mseed[, core.rxn := id %in% core.rxns]
  if(core.only==T)
    mseed <- mseed[core.rxn==T]
  # Remove all duplicate reactions from data base. Keep reactions that are in the core reaction list whereever possible (core.rxns).
  mseed[, rxn.hash := generate_rxn_stoich_hash(stoichiometry, reversibility)]
  mseed <- mseed[order(rxn.hash,-core.rxn)]
  dupl.rxns <- mseed[duplicated(rxn.hash),id]
  mseed <- mseed[!duplicated(rxn.hash)]
  mseed <- mseed[order(id)]
  
  # If selected then only consider reactions with have sequence evidence
  if( core.only ){
    idx <- which( !gsub("_.0$","",mod@react_id) %in% core.rxns & !mod@react_id %in% mod.orig@react_id )
    mod <- rmReact(mod, react=idx)
    #idx2<- which( !mseed$id  %in% core.rxns & !mseed$id %in% gsub("_.0$","",mod.orig@react_id) )
    #mseed <- mseed[-c(idx2),]
  }
  
  # # Delete duplicate dummy reactions in full model
  # cat("Deleting redundant reactions...\n")
  # for(rd in dupl.rxns) {
  #   if(!any(grepl(rd,pres.rxns,fixed = T))) {
  #     mod <- rmReact(mod, react = paste0(rd,"_c0"), rm_met = F)
  #   }
  # }
  
  # Blocking redundant dummy reactions
  red.inds <- which(mod@react_id %in% paste0(dupl.rxns,"_c0"))
  mod@lowbnd[red.inds] <- 0
  mod@uppbnd[red.inds] <- 0
  
  # Create a list of all 'underground' or 'dummy' reactions. = Reactions that are in the reactions database but not (yet)
  # part of the model
  dummy.rxns      <- paste0(mseed[,id],"_c0")
  dummy.rxns.inds <- which(mod@react_id %in% dummy.rxns) # CAUTION: Order of dummy.rxns & dummy.rxns.inds not neccessarily the same. use match instead!
  costly.rxns     <- dummy.rxns[!(dummy.rxns %in% paste0(core.rxns,"_c0"))]
  costly.rxns.inds<- which(mod@react_id %in% costly.rxns) # CAUTION: Order of dummy.rxns & dummy.rxns.inds not neccessarily the same. use match instead!
  
  # change bounds 
  mod@lowbnd[dummy.rxns.inds] <- ifelse(mod@lowbnd[dummy.rxns.inds]==0,0,-dummy.bnd)
  mod@uppbnd[dummy.rxns.inds] <- ifelse(mod@uppbnd[dummy.rxns.inds]==0,0,dummy.bnd)

  # get initial growth
  sol <- optimizeProb(mod.orig)
  gr.orig <- sol@lp_obj
  if(gr.orig >= min.gr*diet.scale){
    warning("Original model is already able to produce the target compound. Nothing to do...")
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                core.rxns = core.rxns,
                growth.rate = 0))
  }
  if(sol@lp_stat != ok){
    warning("FBA on original model did not end successfully! Zero solution is feasibile! But gapfilling will start anyway.")
  }
    
  sol <- optimizeProb(mod)
  gr.dummy <- sol@lp_obj
  if(sol@lp_stat!=ok){
    warning(paste0("Full model is already not able to form. There's no way to successful gap-filling."))
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                core.rxns = core.rxns,
                growth.rate = 0))
  }
  
  #dummy.dt <- data.table()

  # while start here
  
  c.coef <- rep(0,mod@react_num)
  c.coef[dummy.rxns.inds] <- core.weight
  # make non sequence-supported reaction more expensive. This ensures that in the case a core reaction fulfils the same function
  # as a costly reacton, the core reaction will carry the flux
  c.coef[costly.rxns.inds] <- dummy.weight
  ind.ex <- grep("^EX_",mod@react_id)
  c.coef[ind.ex] <- 0
  mtf.scale <- 1
  
  # modj_warm <- sysBiolAlg(mod,
  #                         algorithm = "mtf2",
  #                         costcoeffw = c.coef,
  #                         pFBAcoeff = 1e-6)
  if( sybil::SYBIL_SETTINGS("SOLVER") == "glpkAPI" ) {
    #modj_warm <- sysBiolAlg(mod, algorithm = "mtf", costcoeffw = c.coef)
    modj_warm <- sysBiolAlg(mod,
                            algorithm = "mtf2",
                            costcoeffw = c.coef,
                            pFBAcoeff = 1e-3)
  }
  else {
    #modj_warm <- sysBiolAlg(mod, algorithm = "mtf", costcoeffw = c.coef, scaling = mtf.scale)
    modj_warm <- sysBiolAlg(mod,
                            algorithm = "mtf2",
                            costcoeffw = c.coef,
                            pFBAcoeff = 1e-3,
                            scaling = 1)
  }
  sol.fba <- optimizeProb(modj_warm)
  
  if(sol.fba$stat!=ok){
    warning("pFBA did not end successfully.")
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                core.rxns = core.rxns,
                growth.rate = 0))
  }

  # Retrieve list of utilized dummy reactions and add them to the original/draft model
  ko.dt <- data.table(dummy.rxn = dummy.rxns,
                      d.rxn.ind = dummy.rxns.inds,
                      flux = sol.fba$fluxes[dummy.rxns.inds])
  
  ko.dt <- ko.dt[abs(flux) > 0]
  ko.dt[, core := gsub("_.0","",dummy.rxn) %in% core.rxns]
  cat("Utilized candidate reactions: ",nrow(ko.dt))
  
  if( nrow(ko.dt) == 0){ # no dummy reactions is needed
    warning("No dummy reactions utilized in full model. Nothing to add.")
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                core.rxns = core.rxns,
                growth.rate = 0))
  }
  
  rel.rxns <- gsub("_.0","",ko.dt[, dummy.rxn])
  for(j in rel.rxns) {
    i <- which(mseed$id==j)
    mets  <- unlist(str_split(string = mseed[i,compound_ids],pattern = ";"))
    rxn.info <- str_split(unlist(str_split(string = mseed[i,stoichiometry],pattern = ";")), pattern = ":", simplify = T)
    
    met.comp   <- rxn.info[,3]
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
    
    ind.new.mets <- which(met.ids %in% mod.orig@met_id)
    ind.old.mets <- which(mod.orig@met_id %in% met.ids[ind.new.mets])
    
    met.name[ind.new.mets] <- mod.orig@met_name[ind.old.mets]
    
    mod.orig <- addReact(model = mod.orig, 
                         id = paste0(mseed[i,id],"_c0"), 
                         met = met.ids,
                         Scoef = met.scoef,
                         reversible = is.rev, 
                         metComp = as.integer(met.comp)+1,
                         ub = ifelse(only.backwards, 0, 1000),
                         lb = ifelse(is.rev, -1000, 0),
                         reactName = mseed[i, name], 
                         metName = met.name)
  }
  mod.orig <- add_missing_exchanges(mod.orig)
  
  sol <- optimizeProb(mod.orig)
  #if(sol@lp_stat!=ok | sol@lp_obj < min.obj.val*diet.scale){
  if(sol@lp_stat!=ok | sol@lp_obj < 1e-7){
    warning(paste0("Final model cannot produce enough target even when all candidate reactions are added! obj=", sol@lp_obj, " lp_stat=",sol@lp_stat))
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                core.rxns = core.rxns,
                growth.rate = 0))
  }
  # while stop here
  
  ko.dt[, d.rxn.ind := NULL]
  ko.dt[, rxn.ind := match(dummy.rxn, mod.orig@react_id)]
  ko.dt[, ko.obj := NA_real_]
  ko.dt <- ko.dt[order(core, rnorm(.N))] # Randomization here
  ko.dt[, keep := F]
  
  obj.val <- gr.orig
  # Get reaction essentiality & successively remove reactions:
  # successively remove non-essential gap-fill reactions. Starting with non-core reactions.
  
  dummy.rxn.rm <- c()
  for(i in 1:nrow(ko.dt)) {
    bu.lb <- mod.orig@lowbnd[ko.dt[i,rxn.ind]]
    bu.ub <- mod.orig@uppbnd[ko.dt[i,rxn.ind]]
    
    mod.orig@lowbnd[ko.dt[i,rxn.ind]] <- 0
    mod.orig@uppbnd[ko.dt[i,rxn.ind]] <- 0
    
    sol <- optimizeProb(mod.orig) # feasibiliy already checked (zero solution is possible)
    ko.dt[i, ko.obj := sol@lp_obj]
    obj.val <- sol@lp_obj
    
    # restore bounds
    if(obj.val < min.obj.val*diet.scale) {
      ko.dt[i, keep := T]
      mod.orig@lowbnd[ko.dt[i,rxn.ind]] <- bu.lb
      mod.orig@uppbnd[ko.dt[i,rxn.ind]] <- bu.ub
    }
    # Remember reactions to be removed from model to be ready for iterative use
    if(obj.val >= min.obj.val*diet.scale)
      dummy.rxn.rm <- c(dummy.rxn.rm, ko.dt[i,dummy.rxn])
  }
  # Remove reaction from mod.orig if reaction is not needed (obj.val >0 min.obj.val*diet.scale)
  if( length(dummy.rxn.rm)>0 )
    mod.orig <- rmReact(mod.orig, react = dummy.rxn.rm)
  
  rxns.added <- gsub("_.0","",ko.dt[keep==T,dummy.rxn])
  sol <- optimizeProb(mod.orig)
  if(sol@lp_stat!=ok){
    stop("Final model cannot grow. Something is terribly wrong!")
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                core.rxns = core.rxns,
                growth.rate = 0))
  }
  obj.val <- sol@lp_obj
  
  # Generate output and summary info
  cat("\nGapfill summary:\n")
  cat("Added reactions:      ",length(rxns.added),"\n")
  cat("Added core reactions: ",sum(rxns.added %in% core.rxns),"\n")
  cat("Final growth rate:    ",obj.val/diet.scale,"\n")
  
  return(list(model = mod.orig,
              rxns.added = rxns.added,
              core.rxns = core.rxns,
              growth.rate = obj.val/diet.scale,
              full.mod = mod))
}

sync_full_mod <- function(draft.mod, full.mod) {
  # 1. add Reactions of draft, which are not part of full model, to full model
  trans.rxns <- draft.mod@react_id[!(draft.mod@react_id %in% full.mod@react_id)]
  rind.d <- which(draft.mod@react_id %in% trans.rxns)
  for(i in rind.d) {
    metind <- which(draft.mod@S[,i]!=0) 
    full.mod <- addReact(model = full.mod,
                         id = draft.mod@react_id[i],
                         met = draft.mod@met_id[metind],
                         Scoef = draft.mod@S[metind,i],
                         reversible = draft.mod@react_rev[i],
                         lb = draft.mod@lowbnd[i],
                         ub = draft.mod@uppbnd[i],
                         obj = draft.mod@obj_coef[i],
                         reactName = draft.mod@react_name[i],
                         metName = draft.mod@met_name[metind],
                         metComp = draft.mod@met_comp[metind]
    )
  }
  
  # 2. Apply the same exchange constraints
  ex.inds <- data.table(id = draft.mod@react_id[grepl("^EX_",draft.mod@react_id)], id.d = NA_integer_, id.f = NA_integer_,
                        lb = NA_real_, ub = NA_real_)
  ex.inds[, id.d := match(id,draft.mod@react_id)]
  ex.inds[, id.f := match(id,full.mod@react_id)]
  ex.inds[, lb := draft.mod@lowbnd[id.d]]
  ex.inds[, ub := draft.mod@uppbnd[id.d]]
  
  full.mod@lowbnd[ex.inds$id.f] <- ex.inds$lb
  full.mod@uppbnd[ex.inds$id.f] <- ex.inds$ub
  
  # 3. Set same objective
  full.mod@obj_coef <- rep(0, full.mod@react_num)
  obj.d <- draft.mod@react_id[draft.mod@obj_coef == 1]
  if(length(obj.d) == 0)
    stop("ERR: No objective reaction defined in draft model.")
  
  full.mod@obj_coef[which(full.mod@react_id %in% obj.d)] <- 1
  
  return(full.mod)
}
