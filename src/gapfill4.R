gapfill4 <- function(mod.orig, mod.full, rxn.weights, min.gr = 0.01, bcore = 50,
                     dummy.weight = 100, script.dir, core.only = FALSE, verbose=verbose, gs.origin = NA, rXg.tab, env = "") {

  presRxnWeight <- 0.00001 # Reaction weight for pFBA for reactions already present in model or exchange/demand reactions
  
  # backup model
  mod.orig.bak <- mod.orig
  
  # Load dummy model
  mod <- mod.full
  
  mod <- sync_full_mod(mod.orig, mod)
  
  # get currently present reactions
  pres.rxns <- mod.orig@react_id
  pres.rxns <- gsub("_.*","",pres.rxns)
  
  # Get all reactions (non-duplicates) that are not yet part of the model
  mseed <- fread(paste0(script.dir, "/../dat/seed_reactions_corrected.tsv"), header=T, stringsAsFactors = F)
  mseed <- mseed[gapseq.status %in% c("approved","corrected")]
  mseed <- mseed[!(id %in% pres.rxns)]
  
  mseed <- merge(mseed, rxn.weights, by.x = "id", by.y = "seed", all.x = T)
  mseed[, core.rxn := !is.na(weight)] # reactions without blast hit are not considered core reactions for gapfill
  mseed[is.na(weight), weight := dummy.weight]
  
  mseed[core.rxn == T & bitscore < bcore, core.rxn := F] # reactions with a blast hit bitscore below "bcore" are not considered core reactions for gapfill
  
  if(core.only==T) {
    mseed <- mseed[core.rxn==T]
  }
  # Remove all duplicate reactions from data base. Keep reactions that are in the core reaction list whereever possible (core.rxns).
  if(core.only & nrow(mseed) == 0) {
    warning("There are no more core reaction that could be added to the model. Nothing to do...")
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                rxn.weights = rxn.weights,
                growth.rate = 0))
  }
  mseed[, rxn.hash := generate_rxn_stoich_hash(stoichiometry, reversibility)]
  mseed <- mseed[order(rxn.hash,-core.rxn)]
  dupl.rxns <- mseed[duplicated(rxn.hash),id]
  mseed <- mseed[!duplicated(rxn.hash)]
  mseed <- mseed[order(id)]
  
  # If selected then only consider reactions with have sequence evidence
  if( core.only ){
    idx <- which( !gsub("_.0$","",mod@react_id) %in% mseed[core.rxn==T,id] & !mod@react_id %in% mod.orig@react_id )
    mod <- rmReact(mod, react=idx)
  }
  
  # Blocking redundant dummy reactions
  
  red.inds <- which(mod@react_id %in% paste0(dupl.rxns,"_c0"))
  mod@lowbnd[red.inds] <- 0
  mod@uppbnd[red.inds] <- 0
  
  # get initial growth
  sol <- fba(mod.orig)
  gr.orig <- sol@obj
  if(!is.na(gr.orig) && gr.orig >= min.gr){
    warning("Original model is already able to produce the target compound. Nothing to do...")
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                rxn.weights = rxn.weights,
                growth.rate = 0))
  }
  if(sol@stat %notin% stat){
    warning("FBA on original model did not end successfully – solution is infeasibile. Gapfilling will start anyway.")
  }
  
  sol <- fba(mod)
  gr.dummy <- sol@obj
  if(sol@stat %notin% stat){
    warning(paste0("Full model is already not able to yield a feasible solution. There's no way to successful gap-filling."))
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                rxn.weights = rxn.weights,
                growth.rate = 0))
  }
  
  c.coef.dt <- data.table(id = gsub("_.0$","",mod@react_id))
  c.coef.dt <- merge(c.coef.dt, mseed[,.(id, weight, core.rxn)], by.x="id", by.y="id", all.x = T, sort = F)
  c.coef.dt[is.na(weight), weight := presRxnWeight]
  c.coef.dt[id %in% pres.rxns, weight := presRxnWeight]
  c.coef.dt[, candRxn := id %notin% pres.rxns & !grepl("^EX|^DM|^bio1$",id)] # potential new data base reactions for the model
  
  if(gs.origin==1) {
    max.iter  <- 15
    n.iter    <- 1
    pFBAcoeff <- 1e-3
    tol_bnd   <- COBRAR_SETTINGS("TOLERANCE"); min_tol_bnd <- 1e-9
    gf.success <- FALSE
    mod@lowbnd[which(mod@obj_coef == 1)] <- min.gr # enforce minimum required growth rate
    while(!gf.success) {
      # max iterations reached
      if(n.iter > max.iter) {
        stop("Initial gap-filling using all reactions was not successful although all available reactions were considered. Please check your gap-filling medium or consider reducing the minimum required growth rate.")
      }
      
      # otherwise, let's try
      COBRAR_SETTINGS("TOLERANCE",tol_bnd)
      sol.fba <- pfbaHeuristic(mod, costcoeffw = c.coef.dt$weight,
                               pFBAcoeff = pFBAcoeff)
      
      n.iter  <- n.iter + 1
      flxbm <- sol.fba@fluxes[which(mod@obj_coef == 1)]
      
      if(sol.fba@stat %in% stat && flxbm > tol_bnd) {
        # pFBA says its feasible, but let's check if solution is actually valid
        mod.reduced <- rmReact(mod, react = mod@react_id[which(sol.fba@fluxes == 0)])
        sol.reduced <- fba(mod.reduced)
        flxbm.reduced <- sol.reduced@fluxes[which(mod.reduced@obj_coef == 1)]
        if(sol.reduced@stat %notin% stat && tol_bnd > min_tol_bnd) {
          tol_bnd <- max(tol_bnd / 2, min_tol_bnd)
          cat("No feasible gap-fill solution found due to numeric issues. Adjusting bound tolerance to",tol_bnd,"\n")
          next
        } else if(sol.reduced@stat %notin% stat && tol_bnd <= min_tol_bnd) {
          pFBAcoeff <- pFBAcoeff / 2
          cat("No feasible gap-fill solution found due to numeric issues. Relaxing pFBA coefficient 'c' to",pFBAcoeff,"\n")
          next
        } else {
          gf.success <- TRUE
          next
        }
      } else if (sol.fba@stat %in% stat && flxbm <= tol_bnd) {
        # pFBA says its feasible, but growth rate is zero. Let's try with relaxed pFBA
        pFBAcoeff <- pFBAcoeff / 2
        cat("pFBA feasible but with growth rate of 0. Relaxing pFBA coefficient 'c' to",pFBAcoeff,"\n")
        next
      } 
      
      # pFBA says its infeasible. Let's try with relaxed pFBA.
      pFBAcoeff <- pFBAcoeff / 2
      cat("No feasible gap-fill solution found. Relaxing pFBA coefficient 'c' to",pFBAcoeff,"\n")
    }
  } else {
    sol.fba <- pfbaHeuristic(mod, costcoeffw = c.coef.dt$weight,
                             pFBAcoeff = 1e-3)
  }

  if(sol.fba@stat %notin% stat){
    warning("pFBA did not end successfully.")
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                rxn.weights = rxn.weights,
                growth.rate = 0))
  }
  
  # Retrieve list of utilized dummy reactions and add them to the original/draft model
  ko.dt <- data.table(dummy.rxn = c.coef.dt[candRxn == TRUE,id],
                      d.rxn.ind = which(c.coef.dt$candRxn == TRUE),
                      dummy.weight = c.coef.dt[candRxn == TRUE,weight],
                      flux = sol.fba@fluxes[which(c.coef.dt$candRxn == TRUE)])
  
  
  ko.dt <- ko.dt[abs(flux) > 0]
  ko.dt[, core := gsub("_.0","",dummy.rxn) %in% mseed[core.rxn==T,id]]
  cat("Utilized candidate reactions: ",nrow(ko.dt))
  
  if(nrow(ko.dt) == 0){ # no dummy reactions are needed
    warning("No dummy reactions utilized in full model. Nothing to add.")
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                rxn.weights = rxn.weights,
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
    
    ind.notnew.mets <- which(met.ids %in% mod.orig@met_id)
    notnew.mets.names <- mod.orig@met_name[met_pos(mod.orig, met.ids[ind.notnew.mets])]
    met.name[ind.notnew.mets] <- notnew.mets.names
    
    # get gene associations if any
    #dtg.tmp <- rXg.tab[seed == mseed[i,id] & rm == F, .(complex,gene)]
    #dtg.tmp <- rXg.tab[seed == mseed[i,id] & rm == F][order(-bitscore)][1, .(complex,gene)]
    dtg.tmp <- data.table()
    if(nrow(dtg.tmp)>0)
      gpr.tmp <- get_gene_logic_string(dtg.tmp$complex, dtg.tmp$gene)
    else
      gpr.tmp <- ""
    
    mod.orig <- addReact(model = mod.orig,
                         id = paste0(mseed[i,id],"_c0"),
                         met = met.ids,
                         Scoef = met.scoef,
                         metComp = mod.orig@mod_compart[as.integer(met.comp)+1],
                         ub = ifelse(only.backwards, 0, 1000),
                         lb = ifelse(is.rev, -1000, 0),
                         reactName = mseed[i, name.x],
                         metName = met.name)
    mod.orig@react_attr[which(mod.orig@react_id == paste0(mseed[i,id],"_c0")),c("gs.origin","seed")] <- data.frame(gs.origin = gs.origin,
                                                                                                                   seed = mseed[i,id],
                                                                                                                   stringsAsFactors = F)
    tmp.id <- mseed[i,id]
    if(tmp.id %in% rxn.weights$seed) {
      c.names <- intersect(colnames(rxn.weights), colnames(mod.orig@react_attr))
      mod.orig@react_attr[which(mod.orig@react_id == paste0(mseed[i,id],"_c0")), c.names] <- rxn.weights[seed == tmp.id, ..c.names]
    }
    
  }
  
  # adust environment if needed
  if(env[1] != "")
    mod.orig <- adjust_model_env(mod.orig, env, script.dir)
  
  mod.orig <- add_missing_exchanges(mod.orig)
  
  sol <- fba(mod.orig)
  if(sol@stat %notin% stat | sol@obj < 1e-7){
    warning(paste0("Final model ", mod.orig@mod_name, " cannot produce enough target even when all candidate reactions are added! obj=", sol@obj, " stat=",sol@stat))
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                rxn.weights = rxn.weights,
                growth.rate = 0))
  }
  
  ko.dt[, d.rxn.ind := NULL]
  ko.dt[, rxn.ind := match(dummy.rxn, gsub("_.*","",mod.orig@react_id))]
  ko.dt[, ko.obj := NA_real_]
  #ko.dt <- ko.dt[order(core, rnorm(.N))] # Randomization here
  #ko.dt <- ko.dt[order(!core, dummy.rxn,decreasing = T)] # No Randomization here
  #ko.dt <- ko.dt[order(core, dummy.rxn,decreasing = F)] # No Randomization here
  #ko.dt <- ko.dt[order(core, abs(flux))] # no randomization. Sorting by abs flux # BEST CHOICE as it preferrably removes reactions that carry only small flux
  ko.dt <- ko.dt[order(core, -dummy.weight)]
  ko.dt[, keep := F]
  #if(core.only ==F)
  #  print(ko.dt)
  
  obj.val <- gr.orig
  # Get reaction essentiality & successively remove reactions:
  # successively remove non-essential gap-fill reactions. Only non-core reactions.
  dummy.rxn.rm <- c()
  for(i in 1:nrow(ko.dt[core==FALSE])) {
    bu.lb <- mod.orig@lowbnd[ko.dt[i,rxn.ind]]
    bu.ub <- mod.orig@uppbnd[ko.dt[i,rxn.ind]]
    
    mod.orig@lowbnd[ko.dt[i,rxn.ind]] <- 0
    mod.orig@uppbnd[ko.dt[i,rxn.ind]] <- 0
    
    sol <- fba(mod.orig) # feasibility already checked (zero solution is possible)
    ko.dt[i, ko.obj := sol@obj]
    obj.val <- sol@obj

    # restore bounds
    if(is.na(obj.val) || obj.val < min.gr) {
      ko.dt[i, keep := T]
      mod.orig@lowbnd[ko.dt[i,rxn.ind]] <- bu.lb
      mod.orig@uppbnd[ko.dt[i,rxn.ind]] <- bu.ub
    }
    # Remember reactions to be removed from model to be ready for iterative use
    if(!is.na(obj.val) && obj.val >= min.gr)
      dummy.rxn.rm <- c(dummy.rxn.rm, ko.dt[i,dummy.rxn])
  }
  # Remove reaction from mod.orig if reaction is not needed (obj.val > min.obj.val)
  if( length(dummy.rxn.rm)>0 )
    mod.orig <- rmReact(mod.orig, react = paste0(dummy.rxn.rm,"_c0"))
  
  rxns.added <- gsub("_.0","",ko.dt[keep==T | core==T,dummy.rxn])
  sol <- fba(mod.orig)
  if(sol@stat %notin% stat){
    stop("Final model cannot grow. Something is terribly wrong!")
    return(list(model = mod.orig.bak,
                rxns.added = c(),
                rxn.weights = rxn.weights,
                growth.rate = 0))
  }
  obj.val <- sol@obj
  
  # Generate output and summary info
  cat("\nGapfill summary:\n")
  cat("Added reactions:      ",length(rxns.added),"\n")
  cat("Added core reactions: ",sum(rxns.added %in% ko.dt[core==T,dummy.rxn]),"\n")
  cat("Final growth rate:    ",obj.val,"\n")
  
  return(list(model = mod.orig,
              rxns.added = rxns.added,
              rxn.weights = rxn.weights,
              growth.rate = obj.val,
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
  full.mod@obj_coef <- rep(0, react_num(full.mod))
  obj.d <- draft.mod@react_id[draft.mod@obj_coef == 1]
  if(length(obj.d) == 0)
    stop("ERR: No objective reaction defined in draft model.")
  
  full.mod@obj_coef[which(full.mod@react_id %in% obj.d)] <- 1
  
  return(full.mod)
}


