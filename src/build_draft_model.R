library(getopt)
library(stringr)

build_draft_model_from_reaction_list <- function(rxns, gram, model.name = "draft_mod", script.dir) {
  require(data.table)
  require(stringr)
  require(sybil)
  
  mseed <- fread(paste0(script.dir, "/../dat/seed_reactions_corrected.tsv"), header=T, stringsAsFactors = F)
  mseed <- mseed[gapseq.status %in% c("approved","corrected")]
  mseed <- mseed[order(id)]
  
  mseed <- mseed[id %in% rxns]
  
  #remove duplicate reactions
  mseed[, rxn.hash := generate_rxn_stoich_hash(stoichiometry, reversibility)]
  mseed <- mseed[order(rxn.hash,id)]
  #dupl.rxns <- mseed[duplicated(rxn.hash),id]
  #mseed <- mseed[!duplicated(rxn.hash)]
  
  mod <- modelorg(name = model.name,id = model.name)
  for(i in (1:nrow(mseed))) {
    cat("\r",i,"/",nrow(mseed))
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
    
    ind.new.mets <- which(met.ids %in% mod@met_id)
    ind.old.mets <- which(mod@met_id %in% met.ids[ind.new.mets])
    
    met.name[ind.new.mets] <- mod@met_name[ind.old.mets]
    
    mod <- addReact(model = mod, 
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
  cat("\n")
  
  # Adding Biomass reaction
  if(gram == "neg")
    dt.bm <- fread(paste0(script.dir, "/../dat/seed_biomass.DT_gramNeg.tsv"))
  if(gram == "pos")
    dt.bm <- fread(paste0(script.dir, "/../dat/seed_biomass.DT_gramPos.tsv"))
  
  mod <- addReact(mod,id = "bio1", met = dt.bm$id, Scoef = dt.bm$stoich, reversible = F, lb = 0, ub = 1000, obj = 1, 
                  reactName = paste0("Biomass reaction ",ifelse(gram=="neg","(gram -)","(gram +)")), metName = dt.bm$name, metComp = dt.bm$comp)
  
  mod <- add_missing_exchanges(mod) 
  
  return(mod)
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

source(paste0(script.dir,"/add_missing_exRxns.R"))
source(paste0(script.dir,"/generate_rxn_stoich_hash.R"))

# get options first
spec <- matrix(c(
  'core.reactions', 'c', 1, "character", "List (space-separated) of reactions to be added to the model.",
  'gram', 'g', 1, "character", "Gram \"pos\" OR \"neg\" ?",
  'model.name', 'n', 2, "character", "Name of draft model network. Default: \"draft_mod\"",
  'output.dir', 'o', 2, "character", "Directory to store results. Default: \".\"",
  'sbml.output', 's', 2, "logical", "Should the gapfilled model be saved as sbml? Default: FALSE"
), ncol = 5, byrow = T)

opt <- getopt(spec)

# Help Screen
if ( is.null(opt$core.reactions) | is.null(opt$gram)){
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# Setting defaults if required
if ( is.null(opt$model.name) ) { opt$model.name = "draft_mod" }
if ( is.null(opt$output.dir) ) { opt$output.dir = "." }
if ( is.null(opt$sbml.output) ) { opt$sbml.output = F }

# Arguments:
core.rxn.file <- opt$core.reactions
model.name    <- opt$model.name
gram          <- opt$gram
output.dir    <- opt$output.dir

# Get list of core-reactions from one or more files (given as comma seperated string by -c)
core.rxns <- c()
for( file in unlist(str_split(core.rxn.file, ",")) ){
  input.rxns <- readLines(file)
  core.rxns <- c(core.rxns, unlist(str_split(input.rxns, " ")))
}
core.rxns <- unique(core.rxns)
cat("\nLoading input reactions:", length(core.rxns), "\n")

# construct draft model
mod <- build_draft_model_from_reaction_list(rxns = core.rxns, gram = gram, model.name = model.name, script.dir = script.dir)

# save draft model
saveRDS(mod,file = paste0(model.name, ".RDS"))