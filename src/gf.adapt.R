library(data.table)
library(tools)
library(stringr)

# get current script path
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  # RStudio specific code
  script.dir    <- dirname(rstudioapi::getSourceEditorContext()$path)
} else{
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  script.dir  <- dirname(script.name)
}

# select solver
if( "cplexAPI" %in% rownames(installed.packages()) ){
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
}else{
  warning("glpkAPI is used but cplexAPI is recommended because it is much faster")
  sybil::SYBIL_SETTINGS("SOLVER","glpkAPI"); ok <- 5
}

#setwd(script.dir)

# Arguments:
media.file       <- "../dat/media/TSBmed.csv"
media.org <- fread(paste0(script.dir,"/../dat/media/MM_glu.csv")) # use minimal medium
#media.org <- fread(paste0(script.dir,"/../dat/media/MM_anaerobic_CO2_H2.csv")) # use minimal medium
#media.org <- fread(paste0(script.dir,"/../dat/media/Mineral_salt.csv")) # use minimal medium
#fullmod.file     <- paste0(script.dir, "/../dat/full.model.RDS")
target.met       <- "cpd11416"
bcore            <- 50

# Parameters:
dummy.weight <- 100
min.obj.val  <- 0.05
sbml.export  <- FALSE 

# Get list of core-reactions from one or more files (given as comma seperated string by -c)
#rxn.weights <- readRDS(rxn.weights.file)
#rXg.tab     <- readRDS(rxnXgene.table)
#rxn.blast   <- fread(rxn.blast.file,fill = T)

# Little helpers
source(paste0(script.dir,"/add_missing_exRxns.R"))
source(paste0(script.dir,"/constrain.model.R"))
source(paste0(script.dir,"/gapfill4.R"))
source(paste0(script.dir,"/generate_rxn_stoich_hash.R"))
source(paste0(script.dir,"/get_gene_logic_string.R"))
rm.na <- function(vec){
  idx <- which(is.na(vec))
  if (length(idx) > 0) 
    return(vec[-idx])
  else return(vec)
}

# database files
carbon.source <- fread(paste0(script.dir, "/../dat/subex.tbl"))
meta.db       <- fread(paste0(script.dir, "/../dat/meta_pwy.tbl"))
seed.db       <- fread(paste0(script.dir, "/../dat/seed_reactions_corrected.tsv"))
seed.db.met   <- fread(paste0(script.dir, "/../dat/seed_metabolites_edited.tsv"))


check_lethal <- function(model, rxn_list, med.file=NA, med=NA){
  if(is.na(med.file) & all(is.na(med))) stop("Provide media or media file")
  lethal.dt <- data.table()
  #lethal_list <- sapply(seq_along(rxn_list), function(i){
  for(i in seq_along(rxn_list)){
    if( str_length(rxn_list[i]) > 0){
      rxn.del <- unlist(str_split(rxn_list[i], " "))
      rxn.idx <- match(rxn.del, model@react_id)
      rxn.idx.na <- which(is.na(rxn.idx))
      rxn.idx[rxn.idx.na] <- match(rxn.del[rxn.idx.na], str_extract(model@react_id, "rxn[0-9]+"))
      
      #if( all(is.na(str_extract(rxn.del, "_[a-z][0-9]$"))) ){
      #  rxn.idx <- na.omit(match(rxn.del, str_extract(model@react_id, "rxn[0-9]+")))  
      #}else{
      #  rxn.idx <- na.omit(match(rxn.del, model@react_id))  
      #}
      if( length(na.omit(rxn.idx)) > 0 ){
        if( !is.na(med.file)) mod.tmp <- constrain.model(model, media.file = med.file)
        else mod.tmp <- constrain.model(model, media = med)
        mod.tmp <- rmReact(mod.tmp, react=model@react_id[na.omit(rxn.idx)])
        
        sol <- optimizeProb(mod.tmp, retOptSol=F)
        lethal <- !(sol$stat == ok & sol$obj >= 1e-7)  
        lethal.dt <- rbind(lethal.dt, data.table(nr=i, rxn=paste0(model@react_id[na.omit(rxn.idx)],collapse=","), lethal=lethal))
      } else lethal.dt <- rbind(lethal.dt, data.table(nr=i, rxn="", lethal=NA))
    }else lethal.dt <- rbind(lethal.dt, data.table(nr=i, rxn="", lethal=NA))
  }
  return(lethal.dt)
}

# prepare model for biolog-like test
esp.mode <- function(model, media){
  model <- constrain.model(model, media = media) # constrain model
  model@obj_coef <- rep(0,model@react_num)
  
  # add biolog like test
  mql <- "cpd15499[c0]"; mqn <- "cpd15500[c0]"
  uql <- "cpd15561[c0]"; uqn <- "cpd15560[c0]"
  h   <- "cpd00067[c0]"
  nad  <- "cpd00003[c0]"; nadh  <- "cpd00004[c0]"
  fdox <- "cpd11621[c0]"; fdred <- "cpd11620[c0]" # ferredoxin
  pql  <- "cpd27796[c0]"; pqn   <- "cpd27797[c0]" # plastoquinone
  model <- addReact(model, "ESP1", met=c(mql,h,mqn), Scoef=c(-1,2,1), lb=0, ub=1000, metComp = rep(1,3))
  model <- addReact(model, "ESP2", met=c(uql,h,uqn), Scoef=c(-1,2,1), lb=0, ub=1000, metComp = rep(1,3))
  model <- addReact(model, "ESP3", met=c(nadh,h,nad), Scoef=c(-1,1,1), lb=0, ub=1000, metComp = rep(1,3))                                                         
  model <- addReact(model, "ESP4", met=c(fdred,fdox), Scoef=c(-1,1), lb=0, ub=1000, metComp = rep(1,2))
  model <- addReact(model, "ESP5", met=c(pql,h,pqn), Scoef=c(-1,2,1), lb=0, ub=1000, metComp = rep(1,3))
  model <- changeObjFunc(model, react=c("ESP1", "ESP2", "ESP3", "ESP4", "ESP5"), obj_coef=c(1,1,1,1,1))
  
  return(invisible(model))
}

add_growth <- function(model.orig, add.met.id=NA, weights=NA, genes=NA, verbose=F, only.core=F, fullmod){
  # This here is needed if another draft than GapSeq's own draft networks are gapfilled
  if((!"gs.origin" %in% colnames(model.orig@react_attr))) {
    model.orig@react_attr <- data.frame(seed      = gsub("_.0","",model.orig@react_id),
                                        gs.origin = 0,
                                        stringsAsFactors = F)
  }
  if( !is.na(weights) ) rxn.weights <- readRDS(weights)
  if( !is.na(genes) )   rXg.tab     <- readRDS(genes)
  
  if( nrow(seed.db.met[id==add.met.id]) == 0 ) stop(paste("Metabolite not found in database "),add.met.id)
  add.met.name <- seed.db.met[id==add.met.id, name]
  add.met.ex   <- paste0("EX_",add.met.id, "_e0")
  print(add.met.ex)
  
  if( !add.met.ex %in% model.orig@react_id ){
    print(paste0("Added exchange reaction for:", add.met.id))
    mod.adapt  <- sybil::addReact(model = model.orig,
                                  id = add.met.ex, 
                                  met = paste0(add.met.id, "[e0]"),
                                  Scoef = -1,
                                  reversible = T, 
                                  metComp = 2,
                                  ub = 1000,
                                  lb = 0,
                                  reactName = paste0(add.met.name, " Exchange"), 
                                  metName = add.met.name)
  }else{
    mod.adapt <- model.orig
  }
  
  media <- media.org[name!="D-Glucose"]
  #media <- media.org[!name %in% c("Benzoate", "co2")]
  media <- rbind(media, data.table(compounds=add.met.id, name=add.met.name, maxFlux=100))

  # add metabolite objective + sink
  mod.adapt    <- esp.mode(mod.adapt, media)

  sol <- optimizeProb(mod.adapt, retOptSol=F)
  #cat("Check growth: ", sol$stat==ok, "\t", sol$obj, "\n")
  
  if(sol$stat == ok & sol$obj >= 1e-7){
    warning(paste("Model is already growing with", add.met.id, add.met.name))
    if( verbose ){
      cat("solution", sol$obj, "\n")
      check_growth(mod.adapt, printExchanges = T, useNames = T)
      idx <- which(colSums(mod.adapt@S[grep(add.met.id, mod.adapt@met_id),]) != 0 & sol$fluxes != 0)
      printReaction2(mod.adapt, react=idx, useNames = T)
    } 
    return(invisible(model.orig))
  }else{
    cat("\nTry to gapfill", add.met.name, add.met.id, "\n")
    #invisible(capture.output( 
    mod.adapt.lst <- gapfill4(mod.orig = mod.adapt, 
                              mod.full = fullmod, 
                              rxn.weights = copy(rxn.weights),  
                              min.gr = min.obj.val,
                              bcore = bcore,
                              dummy.weight = dummy.weight,
                              script.dir = script.dir,
                              core.only = only.core,
                              verbose=verbose,
                              gs.origin = 10,
                              rXg.tab = rXg.tab) 
    #))
    new.reactions <- mod.adapt.lst$rxns.added
    if( length(new.reactions) > 0 ){
      #if( verbose ) cat("Added reactions:", new.reactions, "\n")
      mod.adapt <- mod.adapt.lst$model
    }
  }
  #print(creatioExNihilo(mod.adapt.lst$model,short = T))
  #check_growth(mod.adapt.lst$model, medium = c(media$compounds), printExchanges = T, useNames = T, mtf = T)
  
  sol <- optimizeProb(mod.adapt, retOptSol=F)
  #cat(sol$stat==ok, "\t", sol$obj, "\n")
  cat("\n")
  
  check <- check_lethal(mod.adapt, rxn_list = setdiff(mod.adapt@react_id, model.orig@react_id), med = media)
  check$status <- rxn.weights$status[match(str_extract(check$rxn, "rxn[0-9]+"), rxn.weights$seed)]
  print(check)
  mod.adapt <- rmReact(mod.adapt, react=check[lethal==F, rxn])
  rm.idx <- na.omit(match(c("ESP1","ESP2"), mod.adapt@react_id))
  if( length(rm.idx) > 0 ) mod.adapt <- rmReact(mod.adapt, react=rm.idx) 
  mod.adapt <- changeObjFunc(mod.adapt, react=model.orig@react_id[model.orig@obj_coef != 0])
  ex1 <- findExchReact(model.orig); ex2 <- findExchReact(mod.adapt)
  mod.adapt@lowbnd[ex2@react_pos] <- 0
  mod.adapt <- changeBounds(mod.adapt, react = ex1@react_id, lb = ex1@lowbnd) # reset medium 
  rxn.added <- setdiff(mod.adapt@react_id, model.orig@react_id)
  cat("Added reactions:", rxn.added, "\n")
  printReaction(mod.adapt, react=rxn.added)
  return(invisible( mod.adapt ))
}

rm_growth <- function(model.orig, del.met.id, use.media=NA, rxn.blast.file, verbose=F){
  # This here is needed if another draft than GapSeq's own draft networks are gapfilled
  rxn.blast   <- fread(rxn.blast.file, header=T, stringsAsFactors = F, blank.lines.skip = T)
  if((!"gs.origin" %in% colnames(model.orig@react_attr))) {
    model.orig@react_attr <- data.frame(seed      = gsub("_.0","",model.orig@react_id),
                                      gs.origin = 0,
                                      stringsAsFactors = F)
  }
  
  if( nrow(seed.db.met[id==del.met.id]) == 0 ) stop(paste("Metabolite not found in database "),del.met.id)
  del.met.name <- seed.db.met[id==del.met.id, name]
  del.met.ex   <- paste0("EX_",del.met.id, "_e0")
  
  if( !del.met.ex %in% model.orig@react_id ){
    warning(paste0("Model does not have exchange reaction:", del.met.id))
    return(invisible(model.orig))
  }
  if ( !all(is.na(use.media)) ){
    media.user <- use.media
  }else{
    ex <- findExchReact(model.orig); ex <- ex[which(ex@lowbnd<0)]
    media.user <- data.table(compounds=str_extract(ex@met_id, "cpd[0-9]+"),
                             name=ex@react_id,
                             maxFlux=-ex@lowbnd)
  }

  media <- media.org[!name %in% c("D-Glucose", "Glucose")]
  #media <- media.org[!name %in% c("Benzoate", "co2")]
  media <- rbind(media, data.table(compounds=del.met.id, name=del.met.name, maxFlux=100))
  
  model <- esp.mode(model.orig, media)
  sol <- optimizeProb(model, retOptSol=F)
  #cat(sol$stat==ok, "\t", sol$obj, "\n")
  if( !(sol$stat==ok & sol$obj >= 1e-7 ) ){
    warning(paste("Model already cannot grow with:", del.met.id, del.met.name))
    return(invisible(model.orig))
  }else{
    cat("\nTry to remove growth for", del.met.name, del.met.id, "\n")
  }

  del.pat <- paste0(del.met.name, ".*degradation|degradation.*", del.met.name)
  pwy.involved <- meta.db[ (grepl(del.met.name,meta.db$name,ignore.case = T) | grepl(del.met.name,meta.db$altname,ignore.case = T) | grepl(del.met.name,meta.db$hierarchy,ignore.case = T)) & grepl("degradation",meta.db$name,ignore.case = T),id]
  rxn.involved.meta <- unique(unlist(str_split(meta.db[id %in% pwy.involved, reaId], ",")))
  
  idx <- na.omit(match(rxn.involved.meta, rxn.blast$rxn))
  rxn.involved <- rxn.blast[unique(idx),.(rxn,name,bitscore,dbhit)]
  rxn.transporter <- na.omit(str_extract(model@react_id[which(S(model)[match(paste0(del.met.id, "[e0]"), model@met_id),]!=0)], "rxn[0-9]+")) # ATTENTION adding transporter to candidate list!!
  rxn.involved <- rbind(rxn.involved, data.table(rxn="TR", name="transporter", bitscore=NA, dbhit=paste(rxn.transporter, collapse = " ")))
  cat("\tCandidate reactions found:", nrow(rxn.involved), "\n")
  
  #rxn.lethal      <- check_lethal(model, rxn.involved[,dbhit], med=media.user) 
  rxn.lethal      <- check_lethal(model.orig, rxn.involved[,dbhit], med=media.user) # ???
  rxn.stopsgrowth <- check_lethal(model, rxn.involved[,dbhit], med=media)
  rxn.involved$lethal <- rxn.lethal$lethal
  rxn.involved$stopsgrowth <- rxn.stopsgrowth$lethal
  rxn.involved$inmodel <- rxn.stopsgrowth$rxn
  if( verbose ) print(rxn.involved)
  
  # try to find reactions which can be removed to stop growth with substance but which is not essential (i.e. growth with other medium still possilbe)
  rm.candidates <- rxn.involved[stopsgrowth==T & lethal==F]
  cat("\tCandidates for removal:",nrow(rm.candidates), "\n")
  
  # try to take active reactions when nothing else is found @TODO
  if( nrow(rm.candidates) == 0 ){
    warning(paste0("No candidate reactions for deletion found in model"))
    #print(rxn.involved)
    return(invisible(model.orig))
    
    #fva         <- fluxVar(model, percentage=10) 
    #fva.max      <- maxSol(fva, "lp_obj")
    #fva.min      <- minSol(fva, "lp_obj")
    #rxn.active  <- react_id(model)[which(round(fva.max, 3) > 0 | round(fva.min,0) < 0)]
    rxn.active <- react_id(model)[which(round(sol$fluxes, 3) != 0)]
    #rxn.active <- findExchReact(model)@react_id
    cat("\tActive reactions found:", length(rxn.active), "\n")    
    
    # check also active reactions?
    # TODO!!!
    rxn.lethal      <- check_lethal(model, rxn.active, med.file=media.file) # ???
    rxn.stopsgrowth <- check_lethal(model, rxn.active, med=media)
    
    #check_lethal(model, "EX_cpd00076_e0", med=media)
    #check_lethal(model, "rxn09658_c0", med=media)
    #check_growth(rmReact(model, react=c("rxn05655_c0","rxn09658_c0")), medium = media$compounds)
    #printReaction(model, react=which(model@S[grep("cpd00076\\[e0\\]", model@met_id),]!=0))
    #check_growth(model, medium = media$compounds[-22])

    both <- merge(rxn.lethal[,.(rxn,lethal)], rxn.stopsgrowth[,.(rxn,lethal)], by="rxn")
    rm.candidates <- both[lethal.x==F & lethal.y==T]
    rm.rea <- rm.candidates[,rxn]
  }else{
    setorder(rm.candidates, "bitscore")
    #print(head(rm.candidates[,.(rxn, name, inmodel, bitscore, lethal, stopsgrowth)]))
    rm.rea <- unlist(str_split(rm.candidates[1,inmodel],","))
  }
  
  rm.mod <- rmReact(model.orig, react=rm.rea)
  cat("\tDeleted reactions:", paste0(rm.rea, collapse = ", "), "\n")
  #print(optimizeProb(rm.mod))
  
  return(invisible( rm.mod ))
}

check_theo_growth <- function(check.met.id, fullmod){
  if( nrow(seed.db.met[id==check.met.id]) == 0 ) stop(paste("Metabolite not found in database "),check.met.id)
  check.met.name <- seed.db.met[id==check.met.id, name]
  check.met.ex   <- paste0("EX_",check.met.id, "_e0")
  
  media <- media.org[name!="D-Glucose"]
  media <- rbind(media, data.table(compounds=check.met.id, name=check.met.name, maxFlux=100))
  
  # add metabolite objective + sink
  mod.adapt    <- esp.mode(fullmod, media)
  mod.adapt    <- constrain.model(mod.adapt, media=media)
  
  sol <- optimizeProb(mod.adapt, retOptSol=F, algorithm="mtf")
  #sol <- optimizeProb(mod.adapt, retOptSol=F)
  cat("Check growth: ", sol$stat==ok, "\t", sol$obj, "\n")
  
  if(sol$stat == ok & sol$obj >= 1e-7){
    print("Growth possible")
  }else{
    print("No growth possible")
  }
}

addReactSrc <- function(mod, src, rxn.id){
  rxn.idx   <- match(rxn.id, src@react_id)
  met.idx   <- which(src@S[,rxn.idx] != 0)
  met.id    <- src@met_id[met.idx]
  met.scoef <- src@S[met.idx,rxn.idx]
  rxn.rev   <- src@react_rev[rxn.idx]
  rxn.lb    <- src@lowbnd[rxn.idx]
  rxn.ub    <- src@uppbnd[rxn.idx]
  rxn.obj   <- src@obj_coef[rxn.idx]
  
  model.new <- addReact(mod, id=rxn.id, met=met.id,
                       Scoef=met.scoef, reversible=rxn.rev,
                       lb=rxn.lb, ub=rxn.ub, obj=rxn.obj)
  return(model.new)
}


increase_growth <- function(model.orig, min.obj.val, weights=NA, genes=NA, verbose=F, only.core=F, fullmod){
  # This here is needed if another draft than GapSeq's own draft networks are gapfilled
  if((!"gs.origin" %in% colnames(model.orig@react_attr))) {
    model.orig@react_attr <- data.frame(seed      = gsub("_.0","",model.orig@react_id),
                                        gs.origin = 0,
                                        stringsAsFactors = F)
  }
  if( !is.na(weights) ) rxn.weights <- readRDS(weights)
  if( !is.na(genes) )   rXg.tab     <- readRDS(genes)
  
  sol <- optimizeProb(model.orig, retOptSol=F)
  cat("Old growth rate:", round(sol$obj,6), "\n")
  if(sol$stat == ok & sol$obj >= min.obj.val){
    warning(paste("Model already has growth rate >=",min.obj.val))
    return(invisible(model.orig))
  }else{
    cat("\nTry to increase growth", "\n")
    mod.adapt.lst <- gapfill4(mod.orig = model.orig, 
                              mod.full = fullmod, 
                              rxn.weights = copy(rxn.weights),  
                              min.gr = min.obj.val,
                              bcore = bcore,
                              dummy.weight = dummy.weight,
                              script.dir = script.dir,
                              core.only = only.core,
                              verbose=verbose,
                              gs.origin = 10,
                              rXg.tab = rXg.tab) 
  
    new.reactions <- mod.adapt.lst$rxns.added
    if( length(new.reactions) > 0 ){
      #if( verbose ) cat("Added reactions:", new.reactions, "\n")
      mod.adapt <- mod.adapt.lst$model
    }
  }
  sol <- optimizeProb(mod.adapt, retOptSol=F)
  cat("New growth rate:", round(sol$obj,6), "\n")
  rxn.added <- setdiff(mod.adapt@react_id, model.orig@react_id)
  printReaction(mod.adapt, react=rxn.added)
  return(invisible( mod.adapt ))
}
