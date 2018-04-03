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
  'sbml.output', 's', 2, "logical", "Should the gapfilled model be saved as sbml? Default: FALSE"
), ncol = 5, byrow = T)

opt <- getopt(spec)

# Help Screen
if ( !is.null(opt$help) | is.null(opt$model) ){#| length(opt)==1 ) {
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


suppressMessages(library(sybilSBML))
suppressMessages(library(sybil))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(methods))
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1


# Setting defaults if required
if ( is.null(opt$full.model ) ) { opt$full.model = paste0(script.dir,"/dat/full.model.RDS") }
if ( is.null(opt$target.metabolite) ) { opt$target.metabolite = "cpd11416" }
if ( is.null(opt$output.dir) ) { opt$output.dir = "." }
if ( is.null(opt$sbml.output) ) { opt$sbml.output = F }
#if ( is.null(opt$model) ) { opt$model = "" }
if ( is.null(opt$core.reactions) ) { opt$core.reactions = paste0(script.dir,"/dat/myb71-core-Reactions.lst") }
if ( is.null(opt$media) ) { opt$media = paste0(script.dir,"/dat/media/MM_glu.csv") }


# Arguments:
mod.file      <- opt$model
media.file    <- opt$media
fullmod.file  <- opt$full.model
target.met    <- opt$target.metabolite
core.rxn.file <- opt$core.reactions
output.dir    <- opt$output.dir

# Parameters:
diet.scale   <- 4
dummy.bnd    <- 1e-3
dummy.weight <- 100
core.weight  <- 1 
min.obj.val  <- 0.05
sbml.export  <- FALSE 

# Little helpers
source(paste0(script.dir,"/src/add_missing_exRxns.R"))
source(paste0(script.dir,"/src/constrain.model.R"))
source(paste0(script.dir,"/src/gapfill4.R"))
source(paste0(script.dir,"/src/generate_rxn_stoich_hash.R"))
#source("sysBiolAlg_mtfClass2.R")
#source("get_exch_profile.R")

# read full model & target model
mod       <- readRDS(fullmod.file)
mod.orig  <- readSBMLmod(mod.file)
mod.orig  <- add_missing_exchanges(mod.orig)

# constrain model
mod.orig <- constrain.model(mod.orig, media.file = media.file, scaling.fac = diet.scale)
mod.orig@obj_coef <- rep(0,mod.orig@react_num)
#mod      <- constrain.model(mod, media.file = media.file, scaling.fac = diet.scale)

# add metabolite objective + sink
#mod      <- add_met_sink(mod, target.met, obj = 1)
mod.orig <- add_met_sink(mod.orig, target.met, obj = 1)

# Perform gapfill3
cat("\n\n1. Gapfilling with given media using all reactions\n")
mod.fill.lst <- gapfill4(mod.orig = mod.orig, 
                         mod.full = mod, 
                         core.rxn.file = core.rxn.file, 
                         min.gr = min.obj.val,
                         dummy.bnd = dummy.bnd,
                         diet.scale = diet.scale,
                         core.weight = core.weight,
                         dummy.weight = dummy.weight,
                         script.dir = script.dir)

mod.res <- constrain.model(mod.fill.lst$model, media.file = media.file, scaling.fac = 1)
#dt <- get_exch_profile(mod.res)
#dt <- get_active_rxns(mod.res)
#dt[abs(flux)!=0]


if ( TRUE ){
  cat("\n\n2. Gapfilling with minimal medium using only reactions found by sequence\n")

  # constrain model
  mod.orig2 <- mod.res
  media.file2 <- paste0(script.dir,"/dat/media/MM_glu.csv")
  mod.orig2 <- constrain.model(mod.orig2, media.file = media.file2, scaling.fac = diet.scale)
  mod.orig2@obj_coef <- rep(0,mod.orig2@react_num)
  
  bm.ind      <- which(mod.orig2@react_id == "bio1")
  bm.met.inds <- which(mod.orig2@S[,bm.ind]<0)
  bm.met      <- gsub("\\[.0\\]","",mod.orig2@met_id[bm.met.inds])
  bm.met.name <- mod.orig2@met_name[bm.met.inds]
  mod.fill    <- mod.orig2
  
  options(warn=-1)
  for( i in seq_along(bm.met.inds) ){
    target.new <- bm.met[i]
    
    # add metabolite objective + sink
    rm.sink = TRUE
    if( paste0("EX_",target.new,"_c0") %in% react_id(mod.fill) )
      rm.sink = FALSE
    mod.fill  <- add_met_sink(mod.fill, target.new, obj = 1)
    
    sol <- optimizeProb(mod.fill, retOptSol=F)
    
    if(sol$stat == ok & sol$obj >= 1e-6){
      mod.fill@obj_coef <- rep(0,mod.fill@react_num)
    }else{
      #cat("\nTry to gapfill", bm.met.name[i],"\n")
      invisible(capture.output( mod.fill.lst <- gapfill4(mod.orig = mod.fill, 
                                                         mod.full = mod,
                                                         core.rxn.file = core.rxn.file, 
                                                         min.gr = min.obj.val,
                                                         dummy.bnd = dummy.bnd,
                                                         diet.scale = diet.scale,
                                                         core.weight = core.weight,
                                                         dummy.weight = dummy.weight,
                                                         script.dir = script.dir,
                                                         core.only = TRUE) ))
      new.reactions <- mod.fill.lst$rxns.added
      if( length(new.reactions) > 0 ){
        #cat("Added reactions:", new.reactions, "\n")
        mod.fill <- mod.fill.lst$model
      }
      mod.fill@obj_coef <- rep(0,mod.fill@react_num)
    }
    if( rm.sink )
      mod.fill <- rmReact(mod.fill, react=paste0("EX_",target.new,"_c0"))
  }
  options(warn=0)
  
  mod.fill <- changeObjFunc(mod.fill, react=target.met)
  mod.fill <- constrain.model(mod.fill, media.file = media.file, scaling.fac = 1)
  
  cat("Gapfill summary:\n")
  cat("Added reactions:      ",length(mod.fill@react_id)-length(mod.res@react_id),"\n")
  cat("Final growth rate:    ",optimizeProb(mod.fill, retOptSol=F)$obj,"\n")
}






if(!dir.exists(output.dir))
  system(paste0("mkdir ",output.dir))

if(opt$sbml.output)
  writeSBML(mod.fill, filename = paste0(output.dir,"/",gsub(".xml","",basename(mod.file)),"-gapfilled(",target.met,").xml"))

cat(mod.fill.lst$rxns.added, file = paste0(output.dir,"/",gsub(".xml","",basename(mod.file)),"-gapfilled(",target.met,").rxnlst"))

saveRDS(mod.fill, file = paste0(output.dir,"/",gsub(".xml","",basename(mod.file)),"-gapfilled(",target.met,").RDS"))

q(status=0)
