library(getopt)

# get options first
spec <- matrix(c(
  'model', 'm', 1, "character", "GapSeq model to be extended (RDS or SBML)",
  'help' , 'h', 0, "logical", "help",
  'add' , 'a', 2, "character", "reactions or pathways that should be added to the model (comma-separated)",
  'remove' , 'r', 2, "character", "reactions or pathways that should be removed from the model (comma-separated)"
), ncol = 5, byrow = T)

opt <- getopt(spec)

# Help Screen
if ( !is.null(opt$help) | is.null(opt$model)){
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if( is.null(opt$add) & is.null(opt$remove) ){
  print("Please specify reactions or pathway that should be added (-a) or removed (-r)")
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

if( "cplexAPI" %in% installed.packages() )
  suppressMessages(library(cplexAPI))
if( "sybilSBML" %in% installed.packages() )
  suppressMessages(library(sybilSBML))
suppressMessages(library(sybil))
suppressMessages(library(data.table)); setDTthreads(1)
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
#if ( is.null(opt$output.dir) ) { opt$output.dir = "." }
#if ( is.null(opt$sbml.output) ) { opt$sbml.output = F }
if ( is.null(opt$model) ) { opt$model =  paste0(script.dir,"/../toy/myb71.RDS")}

# overwrite -f Option with default (This option might be used in future for a different purpose)

# Arguments:
mod.file            <- opt$model
ids.add             <- opt$add
ids.remove          <- opt$remove

# Parameters:
sbml.export  <- FALSE 

# Little helpers
source(paste0(script.dir,"/add_missing_exRxns.R"))
source(paste0(script.dir,"/addMetAttr.R"))
source(paste0(script.dir,"/addReactAttr.R"))
source(paste0(script.dir,"/construct_full_model.R"))

# database files
meta.pwy <- fread(paste0(script.dir, "/../dat/meta_pwy.tbl"))
custom.pwy <- fread(paste0(script.dir, "/../dat/custom_pwy.tbl"))
all.pwy <- rbind(meta.pwy[,1:9],custom.pwy)
meta.rxn <- fread(paste0(script.dir, "/../dat/meta_rea.tbl"))
seed_x_mets   <- fread(paste0(script.dir,"/../dat/seed_metabolites_edited.tsv"), header=T, stringsAsFactors = F, na.strings = c("null","","NA"))

cat("Loading model files", mod.file, "\n")
mod        <- construct_full_model(script.dir)
if ( toupper(file_ext(mod.file)) == "RDS" ){
  mod.orig <- readRDS(mod.file)
}else{ 
  mod.orig <- readSBMLmod(mod.file)}


#print(ids.add)
#ids <- c("|RIBOSYN2-PWY|", "rxn00023", "14DICHLORBENZDEG-PWY", "rxn05683", "1.1.1.2", "purine hydroxylase", "PYRUVDEH-RXN", "Ammonia-oxidation")

ids2seed <- function(ids){
  id.seed   <- str_extract(ids, "rxn[0-9]+")
  idx.pwy   <- match(gsub("\\|","",ids), gsub("\\|","",all.pwy$id))
  idx.ec    <- str_extract(ids, "[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+")
  idx.rxn   <- match(gsub("\\|","",ids), gsub("\\|","",meta.rxn$id))
  idx.sub   <- unname(sapply(ids, function(id){ hit <- grep(gsub("\\|","",id), all.pwy$hierarchy, ignore.case = T); ifelse(length(hit)>0, paste0(hit, collapse = ","), NA) }))
  id.kegg   <- str_extract(ids, "R[0-9]+")
  
  ids2seed.dt <- data.table()
  for(i in seq_along(ids)){
    if( !is.na(id.seed[i]) ){
      ids2seed.dt <- rbind(ids2seed.dt, data.table(id=ids[i], id.type="seed rxn", db.rxn="", seed=ids[i]))
    }else if ( !is.na(idx.pwy[i]) ){
      rxn <- unlist(str_split(all.pwy[idx.pwy[i], reaId], ","))
      ec  <- unlist(str_split(all.pwy[idx.pwy[i], reaEc], ","))
      rxn.name <- unlist(str_split(all.pwy[idx.pwy[i], reaName], ";"))
      rxn.str <- c()
      for(j in seq_along(rxn)){
        rxn.str <- c(rxn.str, system(paste0(script.dir, "/getDBhit.sh ", paste(rxn[j], paste0("'",rxn.name[j],"'"), ec[j], "seed")), intern=T))
      }
      ids2seed.dt <- rbind(ids2seed.dt, data.table(id=ids[i], id.type="metacyc pwy", db.rxn=rxn, seed=rxn.str))
    } else if ( !is.na(idx.ec[i]) ){
      rxn.str <- system(paste0(script.dir, "/getDBhit.sh ", paste("''", "''", idx.ec[i], "seed")), intern=T)
      ids2seed.dt <- rbind(ids2seed.dt, data.table(id=ids[i], id.type="EC number", db.rxn="", seed=rxn.str))
    }else if ( !is.na(idx.rxn[i]) ){
      rxn <- gsub("\\|","",meta.rxn[idx.rxn[i], id])
      ec  <- meta.rxn[idx.rxn[i], ec]
      rxn.name <- meta.rxn[idx.rxn[i], name]
      rxn.str <- system(paste0(script.dir, "/getDBhit.sh ", paste(rxn, paste0("'",rxn.name,"'"), ec, "seed")), intern=T)
      ids2seed.dt <- rbind(ids2seed.dt, data.table(id=ids[i], id.type="metacyc rxn", db.rxn=ids[i], seed=rxn.str))
    }else if ( !is.na(idx.sub[i]) ){
      idx.sub.split <- as.numeric(unlist(str_split(idx.sub[i], ",")))
      rxn <- unlist(str_split(all.pwy[idx.sub.split, reaId], ","))
      ec  <- unlist(str_split(all.pwy[idx.sub.split, reaEc], ","))
      rxn.name <- unlist(str_split(all.pwy[idx.sub.split, reaName], ";"))
      rxn.str <- c()
      for(j in seq_along(rxn)){
        rxn.str <- c(rxn.str, system(paste0(script.dir, "/getDBhit.sh ", paste(rxn[j], paste0("'",rxn.name[j],"'"), ec[j], "seed")), intern=T))
      }
      ids2seed.dt <- rbind(ids2seed.dt, unique(data.table(id=ids[i], id.type="metacyc sub", db.rxn=rxn, seed=rxn.str)))
    } else if ( !is.na(id.kegg[i]) ){
      rxn.str <- system(paste0(script.dir, "/getDBhit.sh ", paste(id.kegg[i], "''", "''", "seed")), intern=T)
      ids2seed.dt <- rbind(ids2seed.dt, data.table(id=ids[i], id.type="kegg rxn", db.rxn=ids[i], seed=rxn.str))
    }else {
      rxn.str <- system(paste0(script.dir, "/getDBhit.sh ", paste("''", ids[i], "''", "seed")), intern=T)
      ids2seed.dt <- rbind(ids2seed.dt, data.table(id=ids[i], id.type="other", db.rxn="", seed=rxn.str))
    }
  }
  
  return(ids2seed.dt)
}

if( !is.null(ids.add) ){
  rxn.add.dt <- ids2seed(unlist(str_split(ids.add, ",")))
  print(rxn.add.dt)
  rxn.add <- sort(unique(unlist(str_split(rxn.add.dt$seed, " "))))
  rxn.add <- rxn.add[which(rxn.add!="")]
  if( length(rxn.add)==0 ) stop("No model reactions found")
  
  # remove reaction which are already in model
  rxn.new <- setdiff(rxn.add, str_remove(mod.orig@react_id, "_[a-z]0$"))
  print(paste("Reactions already in model: ", paste0(setdiff(rxn.add, rxn.new), collapse = ",")))
  if( length(rxn.new)==0 ) stop("No reactions left.")
  rxn.add <- rxn.new
  
  # remove reacitons which are not in full model
  rxn.avail <- intersect(rxn.add, str_remove(mod@react_id, "_[a-z]0$"))
  print(paste("Reactions not in gapseq reaction database: ", paste0(setdiff(rxn.add, rxn.avail), collapse = ",")))
  if( length(rxn.avail)==0 ) stop("No reactions left.")
  rxn.add <- rxn.avail
  
  mod.out <- add_reaction_from_db(mod.orig, react = rxn.add, gs.origin = 10)
  
  # add metabolite & reaction attributes
  mod.out <- addMetAttr(mod.out, seed_x_mets = seed_x_mets)
  mod.out <- addReactAttr(mod.out)

  print(paste("Added reactions: ", paste0(rxn.add, collapse = ",")))
}

if( !is.null(ids.remove) ){
  rxn.remove.dt <- ids2seed(unlist(str_split(ids.remove, ",")))
  print(rxn.remove.dt)
  rxn.remove <- sort(unique(unlist(str_split(rxn.remove.dt$seed, " "))))
  rxn.remove <- rxn.remove[which(rxn.remove!="")]
  if( length(rxn.remove)==0 ) stop("No model reactions found")
  
  # remove reaction which are already not in model
  rxn.avail  <- intersect(rxn.remove, str_remove(mod.orig@react_id, "_[a-z]0$"))
  if( length(rxn.avail)==0 ) stop("No reactions present in model.")
  rxn.remove <- rxn.avail
  
  if( exists("mod.out") ){
    mod.out <- sybil::rmReact(mod.out, react=rxn.remove)  
  }else{
    mod.out <- sybil::rmReact(mod.orig, react=rxn.remove)  
  }
  print(paste("Removed reactions: ", paste0(rxn.remove, collapse = ",")))
}

out.id <- gsub(".xml|.RDS|.rds","",gsub("-draft","",basename(mod.file)))
out.rds <- paste0("./",out.id,"-adapt",".RDS")
saveRDS(mod.out, file = out.rds)
# Write SBML
source(paste0(script.dir,"/sbml_write.R"))
write_gapseq_sbml(mod.out, paste0(out.id, "-adapt"))
