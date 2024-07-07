construct_full_model <- function(script.dir) {
  suppressMessages(require(data.table))
  suppressMessages(require(stringr))
  suppressMessages(require(cobrar))
  options(warn=1)
  
  # compartment ids and names
  compIDs <- c("c0","e0","p0")
  compNames <- c("Cytosol","Extracellular space","Periplasm")
  
  # Load database 
  db_rxns <- fread(paste0(script.dir, "/../dat/seed_reactions_corrected.tsv"),
                   header=T, stringsAsFactors = F)
  db_rxns <- db_rxns[gapseq.status %in% c("approved","corrected")]
  db_rxns <- db_rxns[order(id)]
  db_mets <- fread(paste0(script.dir,"/../dat/seed_metabolites_edited.tsv"))
  
  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
  # Constuct stoichiometrix matrix S  #
  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
  
  
  # parse stoichiometry
  rxninfos <- lapply(db_rxns$stoichiometry, function(x) {
    tmp <- str_split(unlist(str_split(string = x,pattern = ";")),
                     pattern = ":", simplify = T)
    
    # sometimes a ":" occurs in the metabolite name (e.g. "OPC-8:0), which cause
    # false separation into additional columns. Here: re-concat to single name
    # column
    if(ncol(tmp) > 5) {
      tmp[,5] <- apply(tmp[,5:ncol(tmp)], 1, function(y) paste(y[y != ""],
                                                               collapse = ":"))
      tmp <- tmp[,1:5]
    }
    
    tmp
  })
  
  # add reaction ID
  rxninfos2 <- lapply(1:length(rxninfos), function(i) {
    tmp <- rxninfos[[i]]
    tmp <- cbind(tmp,
                 rep(db_rxns$id[i], nrow(tmp)))
    return(tmp)
  })
  
  # make an easy to handle data.table
  rxninfos <- do.call(rbind, rxninfos2)
  rxninfos <- data.table(rxninfos)
  rm(rxninfos2)
  setnames(rxninfos, c("scoeff","cpd","comp","xyz","cpd.name","rxn.id"))
  
  # correct some Ids and data types
  rxninfos$comp.name <- compIDs[as.integer(rxninfos$comp) + 1L]
  rxninfos[, cpd.id := paste0(cpd,"[",comp.name,"]")]
  rxninfos[, rxn.id := paste0(rxn.id, "_c0")]
  rxninfos[, scoeff := as.double(scoeff)]
  
  cpdinfos <- rxninfos[, .(cpd, cpd.id, cpd.name, comp)]
  cpdinfos <- cpdinfos[!duplicated(cpd.id)]
  cpdinfos[, cpd.name := gsub("^\"|\"$","", cpd.name)]
  cpdinfos[, cpd.name := paste0(cpd.name, "-", compIDs[as.integer(comp) + 1L])]
  tmp_id <- match(cpdinfos$cpd, db_mets$id)
  cpdinfos$charge <- db_mets[tmp_id, charge]
  cpdinfos$chemical.formula <- db_mets[tmp_id, formula]
  
  # table for exchange reactions
  exrxns <- copy(cpdinfos[comp == "1"])
  exrxns[, exid := paste0("EX_", cpd, "_e0")]
  exrxns[, rxnIndex := (1:.N) + nrow(db_rxns)]
  exrxns[, metIndex := match(cpd.id, cpdinfos$cpd.id)]
  exrxns[, exName := paste(cpd.name, "Exchange")]
  
  # construct modelorg object
  mod <- new("ModelOrg", mod_id = "dummy", mod_name = "Full Dummy model with all approved/corrected ModelSEED reactions")
  mod@mod_desc <- "Full Dummy model"
  mod@mod_compart <- compIDs
  mod@mod_compart_name <- compNames
  
  # compounds
  mod@met_id   <- cpdinfos$cpd.id
  mod@met_name <- cpdinfos$cpd.name
  mod@met_comp <- compIDs[as.integer(cpdinfos$comp) + 1L]
  mod@met_attr <- data.frame(charge = cpdinfos$charge,
                             chemicalFormula = cpdinfos$chemical.formula,
                             CVTerms = NA_character_, SBOTerm = NA_character_)
  
  # reactions
  mod@react_id   <- c(paste0(db_rxns$id, "_c0"), exrxns$exid)
  mod@react_name <- c(db_rxns$name, exrxns$exName)
  mod@lowbnd     <- rep(-COBRAR_SETTINGS("MAXIMUM"), react_num(mod))
  mod@uppbnd     <- rep(COBRAR_SETTINGS("MAXIMUM"), react_num(mod))
  mod@obj_coef   <- rep(0, length(mod@react_id))
  mod@uppbnd[db_rxns[,which(reversibility == "<")]] <- 0
  mod@lowbnd[db_rxns[,which(reversibility == ">")]] <- 0
  mod@lowbnd[exrxns$rxnIndex] <- 0
  
  ind_mat_rxns <- cbind(match(rxninfos$cpd.id, mod@met_id),
                        match(rxninfos$rxn.id, mod@react_id),
                        rxninfos$scoeff)
  ind_mat_exch <- cbind(exrxns$metIndex, exrxns$rxnIndex, rep(-1, nrow(exrxns)))
  ind_mat <- rbind(ind_mat_rxns, ind_mat_exch)
  
  # S
  mod@S <- Matrix(0, nrow = met_num(mod), ncol = react_num(mod), sparse = TRUE)
  mod@S[ind_mat[,1:2]] <- ind_mat[,3]
  
  # user constraints
  mod@constraints <- new("Constraints",
                         coeff = as(Matrix(nrow = 0, ncol = ncol(mod@S), sparse = TRUE),
                                    "dMatrix"),
                         lb = numeric(0),
                         ub = numeric(0),
                         rtype = character(0))
  
  # rxnGeneMat adn subsystem matrix
  mod@subSys     <- Matrix(nrow=react_num(mod), ncol = 0, sparse = TRUE)
  
  return(mod)
} 

# TODO: Add ec info and mass balance of reactions in futile cycles
futile_cycle_test <- function(script.dir, env = "") {
  require(cobrar)
  require(data.table)
  require(stringr)
  source(paste0(script.dir, "/print_reaction.R"))
  if("cobrarCPLEX" %in% installed.packages()) {
    library(cobrarCPLEX)
    COBRAR_SETTINGS("SOLVER","cplex")
  }
  
  # parse environment string
  env <- unlist(str_split(env, ","))
  
  # repair path if required
  script.dir <- gsub("\\/$", "", script.dir) # removes tailing "/" if there
  meta.seed <- construct_full_model(script.dir)
  
  
  # adjust environment if needed
  if("highH2" %in% env) {
    env_dt <- fread(paste0(script.dir, "/../dat/env/env_highH2.tsv"), header = F)
    env_dt[, V1 := paste0(V1,"_c0")]
    env_dt <- env_dt[V1 %in% meta.seed@react_id]
    for(i in 1:nrow(env_dt)) {
      if(env_dt[i, V2] == ">")
        meta.seed <- changeBounds(meta.seed, react = env_dt[i, V1],
                                  lb = 0,
                                  ub = COBRAR_SETTINGS("MAXIMUM"))
      if(env_dt[i, V2] == "<")
        meta.seed <- changeBounds(meta.seed, react = env_dt[i, V1],
                                  lb = -COBRAR_SETTINGS("MAXIMUM"),
                                  ub = 0)
      if(env_dt[i, V2] == "=")
        meta.seed <- changeBounds(meta.seed, react = env_dt[i, V1],
                                  lb = -COBRAR_SETTINGS("MAXIMUM"),
                                  ub = COBRAR_SETTINGS("MAXIMUM"))
    }
  }
  
  meta.seed@obj_coef[which(meta.seed@react_id=="rxn00062_c0")] <- 1
  # meta.seed <- changeBounds(meta.seed, "rxn05759_c0", lb = -1000) #  Uncomment to ad-hoc include a futile cycle 
  sol <- pfbaHeuristic(meta.seed)
  
  # save fluxes
  dt <- data.table(id   = meta.seed@react_id, 
                   name = meta.seed@react_name,
                   flux = sol@fluxes
  )
  dt <- dt[abs(flux) > 0.001]
  if(nrow(dt) == 0) {
    return("futile cycle free - Hooray!")
  }
  
  dt[, id.backup := id]
  dt[grepl("^rxn", id), id := gsub("_.0$","", id)]
  
  # get balances for reactions
  #dt[, ec   := ""]
  dt[, CBal := ""]
  #dt[, MBal := ""]
  
  for(i in 1:nrow(dt)) {
    if(grepl("^rxn",dt[i, id])) {
      ind.tmp  <- which(meta.seed@react_id == dt[i, id.backup])
      ind.mets <- which(abs(meta.seed@S[,ind.tmp]) > 0) 
      dt$CBal[i] <- sum(meta.seed@S[ind.mets,ind.tmp] * meta.seed@met_attr$charge[ind.mets])
    }
  }
  
  # Get reaction table 
  gs.rxn <- fread(paste0(script.dir,"/../dat/seed_reactions_corrected.tsv"))
  
  dt <- merge(dt, gs.rxn[, .(id, gapseq.status)], by= "id")
  
  dt$equation <- print_reaction(meta.seed, react = dt$id.backup)
  dt$definition <- print_reaction(meta.seed, react = dt$id.backup,
                                  use.ids = TRUE)
  
  dt[, id.backup := NULL]
  dt[]
  
  fwrite(dt, paste0(script.dir,"/../Futile_cycle.csv"))
  
  return(paste0("At least one futile cycle found. See file: Futile_cycle.csv"))
}
