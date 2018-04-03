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


library(sybilSBML)
library(sybil)
library(data.table)
library(stringr)
library(methods)
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

if(!dir.exists(output.dir))
  system(paste0("mkdir ",output.dir))

if(opt$sbml.output)
  writeSBML(mod.fill.lst$model, filename = paste0(output.dir,"/",gsub(".xml","",basename(mod.file)),"-gapfilled(",target.met,").xml"))

cat(mod.fill.lst$rxns.added, file = paste0(output.dir,"/",gsub(".xml","",basename(mod.file)),"-gapfilled(",target.met,").rxnlst"))


q(status=0)
