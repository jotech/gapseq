library(getopt)

# get options first
spec <- matrix(c(
  'model', 'm', 1, "character", "Model to be gapfilled",
  'help' , 'h', 0, "logical", "help",
  'media', 'n', 1, "character", "tab- or komma separated table for media components. Requires three named columns: 1 - \"compounds\" (for metab. IDs), 2 - \"name\" (metab. name), 3 - \"maxFlux\" (maximum inflow flux)",
  'full.model', 'f', 2, "character","RDS file of the full (dummy) model. (ask Silvio for it :) ). Defaut: full.dummy.model2.RDS",
  'target.metabolite', 't', 2, "character", "ID (without compartment suffix) of metabolite that shall be produced. Default: cpd11416 (Biomass)",
  'core.reactions', 'c', 1, "character", "List (space-separated) of reactions that should preferably be used for gapfilling.",
  'output.dir', 'o', 2, "character", "Directory to store results. Default: \"gapfill\".",
  'sbml.output', 's', 2, "logical", "Should the gapfilled model be saved as sbml? Default: FALSE",
  'verbose', 'v', 2, "logical", "Verbose output and printing of debug messages. Default: FALSE"
), ncol = 5, byrow = T)

opt <- getopt(spec)

# Help Screen
if ( !is.null(opt$help) | is.null(opt$model)){
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# get current script path
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
    # RStudio specific code
    script.dir    <- dirname(rstudioapi::getSourceEditorContext()$path)
} else{
    initial.options <- commandArgs(trailingOnly = FALSE)
    script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
    script.dir  <- dirname(script.name)
}

suppressMessages(library(cplexAPI))
if( "sybilSBML" %in% installed.packages() )
  suppressMessages(library(sybilSBML))
suppressMessages(library(sybil))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(methods))
suppressMessages(library(tools))

# select solver
if( "cplexAPI" %in% rownames(installed.packages()) ){
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
}else{
  warning("glpkAPI is used but cplexAPI is recommended because it is much faster")
  sybil::SYBIL_SETTINGS("SOLVER","glpkAPI"); ok <- 5
}

# Setting defaults if required
if ( is.null(opt$full.model ) ) { opt$full.model = paste0(script.dir,"/dat/full.model.RDS") }
if ( is.null(opt$target.metabolite) ) { opt$target.metabolite = "cpd11416" }
if ( is.null(opt$output.dir) ) { opt$output.dir = "." }
if ( is.null(opt$sbml.output) ) { opt$sbml.output = F }
#if ( is.null(opt$model) ) { opt$model = "" }
if ( is.null(opt$core.reactions) ) { opt$core.reactions = paste0(script.dir,"/dat/myb71-core-Reactions.lst") }
if ( is.null(opt$media) ) { opt$media = paste0(script.dir,"/dat/media/MM_glu.csv") }
if ( is.null(opt$verbose) ) { opt$verbose = F }

# Arguments:
mod.file      <- opt$model
media.file    <- opt$media
fullmod.file  <- opt$full.model
target.met    <- opt$target.metabolite
core.rxn.file <- opt$core.reactions
output.dir    <- opt$output.dir
verbose       <- opt$verbose

# Parameters:
diet.scale   <- 4
dummy.bnd    <- 1e-3
dummy.weight <- 100
core.weight  <- 1 
min.obj.val  <- 0.05
sbml.export  <- FALSE 


# Get list of core-reactions from one or more files (given as comma seperated string by -c)
core.rxns <- c()
for( file in unlist(str_split(core.rxn.file, ",")) ){
  input.rxns <- readLines(file)
  core.rxns <- c(core.rxns, unlist(str_split(input.rxns, " ")))
}
core.rxns <- unique(core.rxns)
cat("\nLoading input reactions:", length(core.rxns), "\n")

# Little helpers
source(paste0(script.dir,"/src/add_missing_exRxns.R"))
source(paste0(script.dir,"/src/constrain.model.R"))
source(paste0(script.dir,"/src/gapfill4.R"))
source(paste0(script.dir,"/src/generate_rxn_stoich_hash.R"))
rm.na <- function(vec){
  idx <- which(is.na(vec))
  if (length(idx) > 0) 
    return(vec[-idx])
  else return(vec)
}


# database files
carbon.source <- fread(paste0(script.dir, "/dat/sub2pwy.csv"))

# read full model & target model
cat("Loading model files", mod.file, "\n")
mod        <- readRDS(fullmod.file)
if ( toupper(file_ext(mod.file)) == "RDS" ){
  mod.orig <- readRDS(mod.file)
}else{ 
  mod.orig <- readSBMLmod(mod.file)}
mod.orig   <- add_missing_exchanges(mod.orig)


# add diffusion reactions
mod.orig       <- add_missing_diffusion(mod.orig)


# create complete medium
if( media.file == "complete" ){
  met.pos <- findExchReact(mod.orig)@met_pos
  met.id  <- gsub("\\[.0\\]","",mod.orig@met_id[met.pos])
  met.name<- mod.orig@met_name[met.pos]
  media <- data.frame(compounds=met.id, name=met.name, maxFlux=100)
  media.file <- paste0(script.dir,"/ALLmed.csv")
  write.csv(media, media.file, quote = F, row.names = F)
}

# constrain model
mod.orig <- constrain.model(mod.orig, media.file = media.file, scaling.fac = diet.scale)
mod.orig@obj_coef <- rep(0,mod.orig@react_num)

# add metabolite objective + sink
mod.orig <- add_met_sink(mod.orig, target.met, obj = 1)

# Perform gapfill
cat("\n\n1. Initial gapfilling: Make model grow on given media using all reactions\n")
mod.fill.lst <- gapfill4(mod.orig = mod.orig, 
                         mod.full = mod, 
                         core.rxn = core.rxn, 
                         min.gr = min.obj.val,
                         dummy.bnd = dummy.bnd,
                         diet.scale = diet.scale,
                         core.weight = core.weight,
                         dummy.weight = dummy.weight,
                         script.dir = script.dir,
                         verbose=verbose)

mod.fill1 <- constrain.model(mod.fill.lst$model, media.file = media.file, scaling.fac = 1)
mod.out <- mod.fill1


if ( TRUE ){
  cat("\n\n2. Biomass gapfilling using core reactions only\n")

  mod.orig2 <- mod.out
  
  # load minimal medium and add available carbon sources
  media2 <- fread(paste0(script.dir,"/dat/media/MM_glu.csv"))
  src.met <- carbon.source[guild %in% c("Carbohydrates", "Polymers", "Carboxylic acids", "Amino acids") & exid_seed %in% mod.orig2@react_id, .(id_seed,name,guild)]
  if( nrow(src.met) == 0)
    stop("No carbon source exchange reactions found in model")
  # if glucose is not usable then add other carbon source(s)
  if( !"alpha-D-Glucose" %in% src.met$name ){
    src.carbo <- src.met[guild=="Carbohydrates"]
    if( nrow(src.carbo)>0 )
      src.add <- src.carbo # if no glucose is there, then add all other available carbohydrates
    else
      src.add <- src.met # if no carbohydrates is avaiable, then take everything else (probably amino acid biosynthesis is not papfilled because amino acids are part of the medium)
    media2 <- rbind(media2, data.table(compounds=gsub("\\[.0\\]","",src.add$id_seed), name=src.add$name, maxFlux=100))  
  }
  
  
  # constrain model  
  mod.orig2 <- constrain.model(mod.orig2, media = media2, scaling.fac = diet.scale)
  mod.orig2@obj_coef <- rep(0,mod.orig2@react_num)
  
  bm.ind      <- which(mod.orig2@react_id == "bio1")
  bm.met.inds <- which(mod.orig2@S[,bm.ind]<0)
  bm.met      <- gsub("\\[.0\\]","",mod.orig2@met_id[bm.met.inds])
  bm.met.name <- mod.orig2@met_name[bm.met.inds]
  mod.fill2    <- mod.orig2
  mod.fill2.counter <- 0
  mod.fill2.names <- c()
  
  if( !verbose ) options(warn=-1)
  for( i in seq_along(bm.met.inds) ){
    cat("\r",i,"/",length(bm.met.inds))
    target.new <- bm.met[i]
    
    # add metabolite objective + sink
    rm.sink = TRUE
    if( paste0("EX_",target.new,"_c0") %in% react_id(mod.fill2) )
      rm.sink = FALSE
    mod.fill2  <- add_met_sink(mod.fill2, target.new, obj = 1)
    
    sol <- optimizeProb(mod.fill2, retOptSol=F)
    
    if(sol$stat == ok & sol$obj >= 1e-6){
      mod.fill2@obj_coef <- rep(0,mod.fill2@react_num)
    }else{
      if( verbose ) cat("\nTry to gapfill", bm.met.name[i],"\n")
      invisible(capture.output( 
                               mod.fill2.lst <- gapfill4(mod.orig = mod.fill2, 
                                                         mod.full = mod,
                                                         core.rxn = core.rxn, 
                                                         min.gr = min.obj.val,
                                                         dummy.bnd = dummy.bnd,
                                                         diet.scale = diet.scale,
                                                         core.weight = core.weight,
                                                         dummy.weight = dummy.weight,
                                                         script.dir = script.dir,
                                                         core.only = TRUE,
                                                         mtf.scale = 2,
                                                         verbose=verbose) ))
      new.reactions <- mod.fill2.lst$rxns.added
      if( length(new.reactions) > 0 ){
        if( verbose ) cat("Added reactions:", new.reactions, "\n")
        mod.fill2 <- mod.fill2.lst$model
        mod.fill2.counter <- mod.fill2.counter + 1
        mod.fill2.names <- c(mod.fill2.names, bm.met.name[i])
      }
      mod.fill2@obj_coef <- rep(0,mod.fill2@react_num)
    }
    if( rm.sink )
      mod.fill2 <- rmReact(mod.fill2, react=paste0("EX_",target.new,"_c0"))
  }
  options(warn=0)
  
  mod.fill2 <- changeObjFunc(mod.fill2, react=paste0("EX_",target.met,"_c0"))
  mod.fill2 <- constrain.model(mod.fill2, media.file = media.file, scaling.fac = 1)
  mod.out <- mod.fill2
  
  cat("\rGapfill summary:\n")
  cat("Filled components:    ",mod.fill2.counter, "(",paste(mod.fill2.names, collapse = ","),")\n")
  cat("Added reactions:      ",length(mod.fill2@react_id)-length(mod.fill1@react_id),"\n")
  cat("Final growth rate:    ",optimizeProb(mod.fill2, retOptSol=F)$obj,"\n")
}


if ( TRUE ){
  cat("\n\n2b. Anaerobic biomass gapfilling using core reactions only\n")
  
  mod.orig2 <- mod.out
  
  # load minimal medium and add available carbon sources
  media2 <- fread(paste0(script.dir,"/dat/media/MM_glu.csv"))
  media2 <- media2[name!="O2"] # remove oxygen
  src.met <- carbon.source[guild %in% c("Carbohydrates", "Polymers", "Carboxylic acids", "Amino acids") & exid_seed %in% mod.orig2@react_id, .(id_seed,name,guild)]
  if( nrow(src.met) == 0)
    stop("No carbon source exchange reactions found in model")
  # if glucose is not usable then add other carbon source(s)
  if( !"alpha-D-Glucose" %in% src.met$name ){
    src.carbo <- src.met[guild=="Carbohydrates"]
    if( nrow(src.carbo)>0 )
      src.add <- src.carbo # if no glucose is there, then add all other available carbohydrates
    else
      src.add <- src.met # if no carbohydrates is avaiable, then take everything else (probably amino acid biosynthesis is not papfilled because amino acids are part of the medium)
    media2 <- rbind(media2, data.table(compounds=gsub("\\[.0\\]","",src.add$id_seed), name=src.add$name, maxFlux=100))  
  }
  
  
  # constrain model  
  mod.orig2 <- constrain.model(mod.orig2, media = media2, scaling.fac = diet.scale)
  mod.orig2@obj_coef <- rep(0,mod.orig2@react_num)
  
  bm.ind      <- which(mod.orig2@react_id == "bio1")
  bm.met.inds <- which(mod.orig2@S[,bm.ind]<0)
  bm.met      <- gsub("\\[.0\\]","",mod.orig2@met_id[bm.met.inds])
  bm.met.name <- mod.orig2@met_name[bm.met.inds]
  mod.fill2    <- mod.orig2
  mod.fill2.counter <- 0
  mod.fill2.names <- c()
  
  if( !verbose ) options(warn=-1)
  for( i in seq_along(bm.met.inds) ){
    cat("\r",i,"/",length(bm.met.inds))
    target.new <- bm.met[i]
    
    # add metabolite objective + sink
    rm.sink = TRUE
    if( paste0("EX_",target.new,"_c0") %in% react_id(mod.fill2) )
      rm.sink = FALSE
    mod.fill2  <- add_met_sink(mod.fill2, target.new, obj = 1)
    
    sol <- optimizeProb(mod.fill2, retOptSol=F)
    
    if(sol$stat == ok & sol$obj >= 1e-6){
      mod.fill2@obj_coef <- rep(0,mod.fill2@react_num)
    }else{
      if( verbose ) cat("\nTry to gapfill", bm.met.name[i],"\n")
      invisible(capture.output( 
        mod.fill2.lst <- gapfill4(mod.orig = mod.fill2, 
                                  mod.full = mod,
                                  core.rxn = core.rxn, 
                                  min.gr = min.obj.val,
                                  dummy.bnd = dummy.bnd,
                                  diet.scale = diet.scale,
                                  core.weight = core.weight,
                                  dummy.weight = dummy.weight,
                                  script.dir = script.dir,
                                  core.only = TRUE,
                                  mtf.scale = 2,
                                  verbose=verbose) ))
      new.reactions <- mod.fill2.lst$rxns.added
      if( length(new.reactions) > 0 ){
        if( verbose ) cat("Added reactions:", new.reactions, "\n")
        mod.fill2 <- mod.fill2.lst$model
        mod.fill2.counter <- mod.fill2.counter + 1
        mod.fill2.names <- c(mod.fill2.names, bm.met.name[i])
      }
      mod.fill2@obj_coef <- rep(0,mod.fill2@react_num)
    }
    if( rm.sink )
      mod.fill2 <- rmReact(mod.fill2, react=paste0("EX_",target.new,"_c0"))
  }
  options(warn=0)
  
  mod.fill2 <- changeObjFunc(mod.fill2, react=paste0("EX_",target.met,"_c0"))
  mod.fill2 <- constrain.model(mod.fill2, media.file = media.file, scaling.fac = 1)
  mod.out <- mod.fill2
  
  cat("\rGapfill summary:\n")
  cat("Filled components:    ",mod.fill2.counter, "(",paste(mod.fill2.names, collapse = ","),")\n")
  cat("Added reactions:      ",length(mod.fill2@react_id)-length(mod.fill1@react_id),"\n")
  cat("Final growth rate:    ",optimizeProb(mod.fill2, retOptSol=F)$obj,"\n")
}




# Add list of exchange reactions for step 3 and 4 in order to check for a wide range of carbon sources or fermentation products
# (Unused exchanges will be deleted afterwards)
idx <- which( !carbon.source$exid_seed %in% mod.out@react_id )
exchanges.new.met  <- rm.na(carbon.source$id_seed[idx])
exchanges.new.name <- rm.na(carbon.source$name[idx])
exchanges.new.ids  <- rm.na(carbon.source$exid_seed[idx])
exchanges.new.used  <- rep(FALSE, length(exchanges.new.ids))  # delete unused addionally added exchange reactions later
mod.out       <- add_exchanges(mod.out, exchanges.new.met, metname=exchanges.new.name)


if ( TRUE ){
  cat("\n\n3. Carbon source gapfilling with core reactions only\n")

  mod.orig3 <- mod.out
  media.org <- fread(paste0(script.dir,"/dat/media/MM_glu.csv")) # use minimal medium
  #media.org <- fread(paste0(script.dir,"/dat/media/Mineral_salt.csv")) # use minimal medium
  
  ex          <- findExchReact(mod.orig3)
  ex.ind      <- ex@react_pos
  ex.id       <- ex@react_id
  ex.met      <- ex@met_id
  ex.met.name <- mod.orig3@met_name[ex@met_pos]
  # Exchange reactions to be ignored (metals etc.)
  ignore <- c("EX_cpd17041_e0", "EX_cpd17042_e0", "EX_cpd17043_e0", "EX_cpd11416_e0", "rxn13782_c0", "rxn13783_c0", "rxn13783_c0", "EX_cpd00001_e0","EX_cpd00007_e0", "EX_cpd00009_e0", "EX_cpd00011_e0" ,"EX_cpd00012_e0", "EX_cpd00030_e0", "EX_cpd00034_e0", "EX_cpd00058_e0", "EX_cpd00063_e0", "EX_cpd00067_e0", "EX_cpd00075_e0","EX_cpd00099_e0", "EX_cpd00149_e0", "EX_cpd00205_e0", "EX_cpd00254_e0", "EX_cpd10515_e0", "EX_cpd00971_e0", "EX_cpd01012_e0", "EX_cpd10516_e0", "EX_cpd11574_e0")
  
  # add metabolite objective + sink
  mod.fill3    <- mod.orig3
  mod.fill3@obj_coef <- rep(0,mod.fill3@react_num)
  
  # add biolog like test
  mql <- "cpd15499[c0]"; mqn <- "cpd15500[c0]"
  uql <- "cpd15561[c0]"; uqn <- "cpd15560[c0]"
  h   <- "cpd00067[c0]"
  #nad <- "cpd00003[c0]"; nadh<- "cpd00004[c0]"
  mod.fill3 <- addReact(mod.fill3, "ESP1", met=c(mql,h,mqn), Scoef=c(-1,2,1), lb=0, ub=1000)
  mod.fill3 <- addReact(mod.fill3, "ESP2", met=c(uql,h,uqn), Scoef=c(-1,2,1), lb=0, ub=1000)
  mod.fill3 <- changeObjFunc(mod.fill3, react=c("ESP1", "ESP2"), obj_coef=c(1,1))
  #mod.fill3 <- addReact(mod.fill3, "ESP3", met=c(nadh,h,nad), Scoef=c(-1,1,1), lb=0, ub=1000)
  #mod.fill3 <- changeObjFunc(mod.fill3, react=c("ESP1", "ESP2", "ESP3"), obj_coef=c(1,1,1))
  mod.fill3.counter <- 0
  mod.fill3.names <- c()
  
  if( !verbose ) options(warn=-1)
  for( i in seq_along(ex.met) ){
    cat("\r",i,"/",length(ex.met))
    if( ex.id[i] %in% ignore ) 
      next
    
    src.met      <- ex.met[i]
    src.met.name <- ex.met.name[i]
    media <- media.org[name!="D-Glucose"]
    #media <- media.org[!name %in% c("Benzoate", "co2")]
    media <- rbind(media, data.table(compounds=gsub("\\[.0\\]","",src.met), name=src.met.name, maxFlux=100))
    
    # constrain model
    mod.fill3 <- constrain.model(mod.fill3, media = media, scaling.fac = diet.scale)
    
    sol <- optimizeProb(mod.fill3, retOptSol=F)
    
    if(sol$stat == ok & sol$obj >= 1e-7){
      #mod.fill3@obj_coef <- rep(0,mod.fill3@react_num)
    }else{
      if( verbose ) cat("\nTry to gapfill", src.met.name, ex@react_id[i], "\n")
      invisible(capture.output( mod.fill3.lst <- gapfill4(mod.orig = mod.fill3, 
                                                         mod.full = mod, 
                                                         core.rxn = core.rxn, 
                                                         min.gr = min.obj.val,
                                                         dummy.bnd = dummy.bnd,
                                                         diet.scale = diet.scale,
                                                         core.weight = core.weight,
                                                         dummy.weight = dummy.weight,
                                                         script.dir = script.dir,
                                                         core.only = TRUE,
                                                         verbose=verbose) ))
      new.reactions <- mod.fill3.lst$rxns.added
      if( length(new.reactions) > 0 ){
        if( verbose ) cat("Added reactions:", new.reactions, "\n")
        mod.fill3 <- mod.fill3.lst$model
        mod.fill3.counter <- mod.fill3.counter + 1
        mod.fill3.names <- c(mod.fill3.names, src.met.name)
        if( ex@react_id[i] %in% exchanges.new.ids) # delete unused addionally added exchange reactions later
          exchanges.new.used[match(ex@react_id[i], exchanges.new.ids)] <- TRUE
      }
    }
  }
  options(warn=0)
  
  mod.fill3 <- rmReact(mod.fill3, react=c("ESP1","ESP2"))
  mod.fill3 <- changeObjFunc(mod.fill3, react=paste0("EX_",target.met,"_c0"))
  mod.fill3 <- constrain.model(mod.fill3, media.file = media.file, scaling.fac = 1)
  mod.out <- mod.fill3
  cat("\rGapfill summary:\n")
  cat("Filled components:    ",mod.fill3.counter, "(",paste(mod.fill3.names, collapse = ","),")\n")
  cat("Added reactions:      ",length(mod.fill3@react_id)-length(mod.fill2@react_id),"\n")
  cat("Final growth rate:    ",optimizeProb(mod.fill3, retOptSol=F)$obj,"\n")
}


if ( TRUE ){
  cat("\n\n4. Checking for fermentation products with core reactions only\n")
  
  mod.orig4 <- mod.out
  media.org <- fread(paste0(script.dir,"/dat/media/MM_glu.csv")) # use minimal medium
  
  ex          <- findExchReact(mod.orig4)
  ex.ind      <- ex@react_pos
  ex.id       <- ex@react_id
  ex.met      <- ex@met_id
  ex.met.name <- mod.orig4@met_name[ex@met_pos]
  # Exchange reactions to be ignored (metals etc.)
  ignore <- c("EX_cpd17041_e0", "EX_cpd17042_e0", "EX_cpd17043_e0", "EX_cpd11416_e0", "rxn13782_c0", "rxn13783_c0", "rxn13783_c0", "EX_cpd00001_e0","EX_cpd00007_e0", "EX_cpd00009_e0", "EX_cpd00011_e0" ,"EX_cpd00012_e0", "EX_cpd00030_e0", "EX_cpd00034_e0", "EX_cpd00058_e0", "EX_cpd00063_e0", "EX_cpd00067_e0", "EX_cpd00075_e0","EX_cpd00099_e0", "EX_cpd00149_e0", "EX_cpd00205_e0", "EX_cpd00254_e0", "EX_cpd10515_e0", "EX_cpd00971_e0", "EX_cpd01012_e0", "EX_cpd10516_e0", "EX_cpd11574_e0")
  
  # add metabolite objective + sink
  mod.fill4    <- mod.orig4
  mod.fill4 <- constrain.model(mod.fill4, media.file = media.file, scaling.fac = diet.scale)
  
  mod.fill4.counter <- 0
  mod.fill4.names <- c()
  
  if( !verbose ) options(warn=-1)
  for( i in seq_along(ex.met) ){
    cat("\r",i,"/",length(ex.met))
    if( ex.id[i] %in% ignore ) 
      next
    
    src.met      <- ex.met[i]
    src.met.name <- ex.met.name[i]
    src.id       <- ex.id[i]
    mod.fill4@obj_coef <- rep(0,mod.fill4@react_num)
    mod.fill4 <- changeObjFunc(mod.fill4, react=src.id, obj_coef=1)
    
    sol <- optimizeProb(mod.fill4, retOptSol=F)
    
    if(sol$stat == ok & sol$obj >= 1e-7){
      #mod.fill4@obj_coef <- rep(0,mod.fill4@react_num)
    }else{
      if( verbose ) cat("\nTry to gapfill", src.met.name, src.id, "\n")
      invisible(capture.output( mod.fill4.lst <- gapfill4(mod.orig = mod.fill4, 
                                                          mod.full = mod, 
                                                          core.rxn = core.rxn, 
                                                          min.gr = min.obj.val,
                                                          dummy.bnd = dummy.bnd,
                                                          diet.scale = diet.scale,
                                                          core.weight = core.weight,
                                                          dummy.weight = dummy.weight,
                                                          script.dir = script.dir,
                                                          core.only = TRUE,
                                                          verbose=verbose) ))
      new.reactions <- mod.fill4.lst$rxns.added
      if( length(new.reactions) > 0 ){
        if( verbose ) cat("Added reactions:", new.reactions, "\n")
        mod.fill4 <- mod.fill4.lst$model
        mod.fill4.counter <- mod.fill4.counter + 1
        mod.fill4.names <- c(mod.fill4.names, src.met.name)
        if( ex@react_id[i] %in% exchanges.new.ids) # delete unused addionally added exchange reactions later
          exchanges.new.used[match(ex@react_id[i], exchanges.new.ids)] <- TRUE
      }
    }
  }
  options(warn=0)
  
  mod.fill4 <- changeObjFunc(mod.fill4, react=paste0("EX_",target.met,"_c0"))
  mod.fill4 <- constrain.model(mod.fill4, media.file = media.file, scaling.fac = 1)
  mod.out <- mod.fill4
  cat("\rGapfill summary:\n")
  cat("Filled components:    ",mod.fill4.counter, "(",paste(mod.fill4.names, collapse = ","),")\n")
  cat("Added reactions:      ",length(mod.fill4@react_id)-length(mod.fill3@react_id),"\n")
  cat("Final growth rate:    ",optimizeProb(mod.fill4, retOptSol=F)$obj,"\n")
}

# delete unused addionally added exchange reactions later
exchanges.rm <- exchanges.new.ids[!exchanges.new.used]
if( length(exchanges.rm) > 0 )
  mod.out <- rmReact(mod.out, react=exchanges.rm)


if(!dir.exists(output.dir))
  system(paste0("mkdir ",output.dir))

if(opt$sbml.output)
  writeSBML(mod.out, filename = paste0(output.dir,"/",gsub(".xml","",basename(mod.file)),"-gapfilled.xml"))

mod.out.rxns.added <- setdiff(mod.out@react_id, mod.orig@react_id)
cat(mod.out.rxns.added, file = paste0(output.dir,"/",gsub(".xml","",basename(mod.file)),"-gapfilled.rxnlst"))

mod.out.rxns.added.without.seq <- setdiff(gsub("_.0","",mod.out.rxns.added), core.rxns)
cat(mod.out.rxns.added.without.seq, file = paste0(output.dir,"/",gsub(".xml","",basename(mod.file)),"-gapfilled.without.seq.rxnlst"))

saveRDS(mod.out, file = paste0(output.dir,"/",gsub(".xml","",basename(mod.file)),"-gapfilled.RDS"))

q(status=0)
