library(getopt)
suppressMessages(library(stringr))
suppressMessages(library(stringi))
suppressMessages(library(R.utils))
library(methods)
suppressMessages(library(IRanges))
options(error=traceback)

# test data
#blast.res <- "../gapseq.publication/gene.essent/GCF_000005845.2_ASM584v2_genomic-all-Reactions.tbl"
#gram <- "auto"
#genome.seq <- "../gapseq.publication/gene.essent/GCF_000005845.2_ASM584v2_genomic.fna.gz"
#high.evi.rxn.BS <- 200
#transporter.res <- "../gapseq.publication/gene.essent/GCF_000005845.2_ASM584v2_genomic-Transporter.tbl"

build_draft_model_from_blast_results <- function(blast.res, transporter.res, biomass = "auto", model.name = NA, genome.seq = NA, 
                                                 script.dir, high.evi.rxn.BS = 200, pathway.pred = NA, min.bs.for.core = 50, 
                                                 curve.alpha = 1) {
  suppressMessages(require(data.table))
  setDTthreads(1)
  suppressMessages(require(stringr))
  #suppressMessages(require(sybil))

  source(paste0(script.dir,"/gram_by_network.R"))
  
  if(is.na(model.name))
    model.name <- gsub("-all-Reactions.tbl","",basename(blast.res), fixed = T)
  
  if(biomass=="auto") {
    if(grepl("\\.gz$", genome.seq)) {
      suppressMessages(require(R.utils))
      genome.seq <- R.utils::gunzip(genome.seq, remove = F, temporary = T, overwrite = T)
    } else {
      
    }

    biomass <- system(paste0(script.dir,"/./predict_biomass_from16S.sh ",genome.seq), intern = T)
    
    cat("\nPredicted biomass: ",biomass,"\n")
    if(biomass == "ambiguous" & !is.na(pathway.pred)) {
        cat("Trying to predict biomass by metabolic network similarity\n")
      biomass <- find_gram_by_network(pathway.pred, script.dir)
        cat("New predicted biomass: ",biomass,"\n")
    }
    
    if(!biomass %in% c("pos","neg", "archaea", "Archaea", "Gram_pos", "Gram_neg")) {
      stop("ERROR: Gram-staining prediction failed or ambiguous result. Please check whether genome sequence contains 16S rRNA gene(s).")

    }
    
  }
  
  seed_x_ec <- fread(paste0(script.dir,"/../dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv"), header=T)
  seed_x_name <- fread(paste0(script.dir,"/../dat/seed_reactions_corrected.tsv"), header=T, stringsAsFactors = F)
  #seed_x_metCyc <- fread(paste0(script.dir,"/../dat/mnxref_seed-other.tsv"), header = T)
  seed_x_aliases <- fread(paste0(script.dir,"/../dat/seed_Enzyme_Name_Reactions_Aliases.tsv"), header=T)
  seed_x_mets   <- fread(paste0(script.dir,"/../dat/seed_metabolites_edited.tsv"), header=T, stringsAsFactors = F, na.strings = c("null","","NA"))
  
  # Prepare candidate reaction tables for draft network and gapfilling
  dt.cand.tmp <- prepare_candidate_reaction_tables(blast.res, transporter.res, high.evi.rxn.BS, min.bs.for.core, curve.alpha)
  dt      <- dt.cand.tmp$dt
  dt.cand <- dt.cand.tmp$dt.cand
  
  dt[!is.na(complex) & is.na(complex.status), complex.status := 0] # incomplete complexes
  
  mseed <- seed_x_name
  mseed <- mseed[gapseq.status %in% c("approved","corrected")]
  mseed <- mseed[order(id)]
  
  # Get all reactions / transporters that have either high sequence evidence (bitscore)
  dt_seed_single_and_there <- copy(dt[(bitscore >= high.evi.rxn.BS & status != "bad_blast") | pathway.status %in% c("full","threshold","keyenzyme")])
  dt_seed_single_and_there <- dt_seed_single_and_there[complex.status != 0 | is.na(complex.status)]
  dt_seed_single_and_there <- dt_seed_single_and_there[order(seed,-bitscore)]
  dt_seed_single_and_there <- dt_seed_single_and_there[!duplicated(seed)]
  
  # check if contig names match conventions and are unique. if not, assign new ones.
  #dt <- fread("Clostridium_difficile_NAP07-all-Reactions.tbl", fill = T)
  contig.names.full  <- unique(dt$stitle)
  contig.names.full  <- contig.names.full[contig.names.full != ""]
  n.contigs          <- length(contig.names.full)
  contig.names.short <- gsub(" .*$","", contig.names.full)
  if(length(unique(contig.names.short)) != n.contigs) {
    warning("Sequence identifiers do not match NCBI standard. Trying to find unique parts in sequence id strings...")
    # try to find a unique parts in contig identifiers
    contig.names.full.tmp <- gsub(":","",contig.names.full, fixed = T)
    id.splits <- str_split(contig.names.full.tmp, pattern = " ")
    nr.splits <- min(unlist(lapply(id.splits, length))) # how many times can we split the identifier string at space (minimum)
    new.id <- c()
    if(nr.splits > 1) {
      for(i in 2:nr.splits) {
        tmp.id <- unlist(lapply(id.splits, function(x) x[i]))
        if(length(unique(tmp.id)) == n.contigs){
          new.id <- tmp.id
          break;
        }
      }
    }
    if(length(new.id) == 0 ) {
      warning("No unique parts in sequence id strings found. Assigning new contig names: \"contig_[i]\"")
      new.id <- paste0("contig_",1:n.contigs)
    }
    dt.new.ids <- data.table(old.contig = contig.names.full, new.contig = new.id)
    for(i in 1:nrow(dt)) {
      cur.stitles <- dt[i,stitle]
      if(cur.stitles != "" & !is.na(cur.stitles)) {
        new.stitle <- dt.new.ids[old.contig == cur.stitles, new.contig]
        dt[i, stitle := new.stitle]
      }
    }
  }
  
  
  # create gene list and attribute table
  cat("Creating Gene-Reaction list... ")

  dt_genes <- copy(dt[!is.na(bitscore)])
  dt_genes[, gene := paste0(gsub(" .*","",stitle),"_",paste(sstart,send, sep = ":"))]

  all_cont <- unique(dt_genes$stitle) # Contigs
  all_cont <- all_cont[all_cont != ""]

  dt_genes[sstart > send, istart := send]
  dt_genes[sstart > send, iend   := sstart]
  dt_genes[sstart < send, istart := sstart]
  dt_genes[sstart < send, iend   := send]
  
  #setkey(dt_genes, stitle, seed)
  
  for(i_cont in all_cont) {
    setkey(dt_genes, stitle, seed)
    dt_genes_cont <- dt_genes[.(i_cont),]
    
    dt_ttt <- dt_genes_cont[,.(gene, istart, iend, bitscore)]
    dt_ttt <- dt_ttt[!duplicated(gene)]
    dt_ttt <- dt_ttt[order(-bitscore)]
    
    ira <- IRanges(start = dt_ttt$istart, end = dt_ttt$iend, names = dt_ttt$gene)
    group <- findOverlaps(ira, ira, type = "within", select = "first")
    dt_ttt$gene2 <- dt_ttt$gene[group]
    
    dt_genes <- merge(dt_genes, dt_ttt[,.(gene, gene2)], by = "gene", all.x = T, 
                      sort = F) 
    dt_genes[!is.na(gene2) & stitle == i_cont, gene := gene2]
    dt_genes[, gene2 := NULL]
  }
  
  dt_genes[, istart := NULL]
  dt_genes[, iend   := NULL]
  
  dt_genes <- dt_genes[!duplicated(paste0(stitle, gene, seed, complex, sep = "$"))]
  
  
  dt_genes[grepl("Subunit undefined", complex), complex := "ZSubunit undefined"]
  dt_genes <- dt_genes[order(stitle,seed,complex,-bitscore)]
  dt_genes <- dt_genes[!duplicated(paste(stitle,seed,gene,sep = "$"))]
  dt_genes[grepl("ZSubunit undefined", complex), complex := "Subunit undefined"]
  
  dt_genes[, rm := F]
  dt_genes$itmp <- 1:nrow(dt_genes)
  setkey(dt_genes, stitle, seed)
  
  for(i_cont in all_cont) {
    all_seed <- unique(dt_genes[.(i_cont), seed]) # Individual reactions on contig
    all_seed <- intersect(mseed$id, all_seed)
    for(i_seed in all_seed) {
      dt_g_tmp <- dt_genes[.(i_cont,i_seed)]
      if(nrow(dt_g_tmp)>1) {
        for(i in 1:(nrow(dt_g_tmp)-1)) {
          astart <- dt_g_tmp[i, sstart]
          aend   <- dt_g_tmp[i, send]
          ind.rm <- dt_g_tmp[(i+1):nrow(dt_g_tmp)][calc_seq_overlap(astart, aend, sstart, send) > 0.5, itmp]
          dt_genes[ind.rm, rm := T]
        }
      }
    }
  }
  cat(length(unique(dt_genes[rm == F & seed %in% mseed$id, paste(gene, sep="$")])),"unique genes on",length(unique(dt_genes[rm == F, stitle])),"genetic element(s)\n")
  
  # create subsys list and attribute table
  dt_subsys <- copy(dt[!is.na(pathway)])
  subsys_unique <- unique(dt_subsys$pathway)
  
  rxns <- unique(dt_seed_single_and_there[,seed])
  
  mseed <- mseed[(id %in% rxns)]
  
  #remove duplicate reactions
  mseed[, rxn.hash := generate_rxn_stoich_hash(stoichiometry, reversibility)]
  mseed <- mseed[order(id, rxn.hash)]
  mseed <- mseed[!duplicated(rxn.hash)]
  
  cat("Constructing draft model: \n")
  mod <- sybil::modelorg(name = model.name, id = model.name)
  mod@mod_desc <- model.name
  mod@react_attr <- data.frame(rxn = character(0), name = character(0), ec = character(0), tc = character(0), qseqid = character(0),
                               pident = numeric(0), evalue = numeric(0), bitscore = numeric(0), qcovs = numeric(0),
                               stitle = character(0), sstart = numeric(0), send = numeric(0), pathway = character(0),
                               status = character(0), pathway.status = character(0), complex = character(0), exception = numeric(0),
                               complex.status = numeric(0), seed = character(0), stringsAsFactors = F)
  mod@subSys <- Matrix::Matrix(F,nrow = 0, ncol = length(subsys_unique),sparse = T)
  colnames(mod@subSys) <- subsys_unique
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
    
    # get reaction-associated stretches of DNA 
    dtg.tmp <- dt_genes[seed == mseed[i,id] & bitscore >= high.evi.rxn.BS & rm == F, .(complex,gene)]
    if(nrow(dtg.tmp)>0)
      gpr.tmp <- get_gene_logic_string(dtg.tmp$complex, dtg.tmp$gene)
    else
      gpr.tmp <- ""
    
    # get reaction-associated subsystems
    dts.tmp <- dt_subsys[seed == mseed[i,id], pathway]
    dts.tmp <- unique(dts.tmp)
    
    mod <- sybil::addReact(model = mod, 
                           id = paste0(mseed[i,id],"_c0"), 
                           met = met.ids,
                           Scoef = met.scoef,
                           reversible = is.rev, 
                           metComp = as.integer(met.comp)+1,
                           ub = ifelse(only.backwards, 0, sybil::SYBIL_SETTINGS("MAXIMUM")),
                           lb = ifelse(is.rev, -sybil::SYBIL_SETTINGS("MAXIMUM"), 0),
                           reactName = mseed[i, name], 
                           metName = met.name,
                           gprAssoc = gpr.tmp,
                           subSys = dts.tmp)
    if(mseed[i,id] %in% dt_seed_single_and_there[,seed])
      mod@react_attr[which(mod@react_id == paste0(mseed[i,id],"_c0")),] <- as.data.frame(dt_seed_single_and_there[seed == mseed[i,id]])
  }
  mod@react_attr$gs.origin <- 0
  mod@react_attr$gs.origin[mod@react_attr$bitscore < high.evi.rxn.BS] <- 9 # Added due to Pathway Topology criteria
  #mod <- add_reaction_from_db(mod, react = c("rxn13782","rxn13783","rxn13784"), gs.origin = 6) # Adding pseudo-reactions for Protein biosynthesis, DNA replication and RNA transcription
  mod <- add_missing_diffusion(mod)
  
  # in case of butyryl CoA:acetate CoA transferase presense - add but transporter as well
  if(any(grepl("rxn00875", mod@react_id)) & all(!grepl("rxn05683", mod@react_id))) {
    mod <- add_reaction_from_db(mod, react = "rxn05683", gs.origin = 1)
  }
  
  #mod@genes_table <- copy(dt[bitscore>0])
  cat("\n")
  warnings()
  # Adding Biomass reaction
  if(biomass == "neg" | biomass == "Gram_neg"){
    ls.bm <- parse_BMjson(paste0(script.dir, "/../dat/biomass/biomass_Gram_neg.json"), seed_x_mets)
    dt.bm <- ls.bm$bmS
    #dt.bm <- fread(paste0(script.dir, "/../dat/biomass/seed_biomass.DT_gramNeg.tsv"))
  }
  if(biomass == "pos" | biomass == "Gram_pos"){
    ls.bm <- parse_BMjson(paste0(script.dir, "/../dat/biomass/biomass_Gram_pos.json"), seed_x_mets)
    dt.bm <- ls.bm$bmS
    #dt.bm <- fread(paste0(script.dir, "/../dat/biomass/seed_biomass.DT_gramPos.tsv"))
  }
  if(biomass == "archaea" | biomass == "Archaea"){
    ls.bm <- parse_BMjson(paste0(script.dir, "/../dat/biomass/biomass_archaea.json"), seed_x_mets)
    dt.bm <- ls.bm$bmS
    #dt.bm <- fread(paste0(script.dir, "/../dat/biomass/seed_biomass.DT_archaea.tsv"))
  }
  
  
  if(biomass %in% c("neg","pos","Gram_neg","Gram_pos")) {
    
    # remove menaquinone8 from anaerobic Biomass if ne novo biosynthesis pathway is absent
    
    if( is.na(dt[grepl("MENAQUINONESYN-PWY",pathway),pathway.status][1]) | is.na(dt[grepl("PWY-5852",pathway),pathway.status][1]) | is.na(dt[grepl("PWY-5837",pathway),pathway.status][1]) ){
      cofac_mass_0 <- dt.bm[met_group == "Cofactors", sum(-Scoef * mass(formula) / 1000)]
      
      dt.bm <- dt.bm[id != "cpd15500[c0]"] # Menaquinone-8
      dt.bm <- dt.bm[id != "cpd15352[c0]"] # 2-Demethylmenaquinone-8
      
      # When removing mets from the biomass reaction we need to rescale the remaining
      # metabolite's coefficients to remain at a netto balance of 1 g/DW.
      # This is done here:
      cofac_mass_1 <- dt.bm[met_group == "Cofactors", sum(-Scoef * mass(formula) / 1000)]
      dt.bm[met_group == "Cofactors", Scoef := Scoef * (cofac_mass_0 / cofac_mass_1)]
    }
  }
  
  mod <- sybil::addReact(mod,id = "bio1", 
                         met = dt.bm$id, 
                         Scoef = dt.bm$Scoef, 
                         reversible = F, 
                         lb = 0, ub = sybil::SYBIL_SETTINGS("MAXIMUM"), 
                         obj = 1, 
                         reactName = ls.bm$name,
                         metName = dt.bm$name,
                         metComp = ifelse(dt.bm$comp=="c",1,ifelse(dt.bm$comp=="e",2,3)))
  mod@react_attr[which(mod@react_id == "bio1"),c("gs.origin","seed")] <- data.frame(gs.origin = 6, seed = "bio1", stringsAsFactors = F)

  # add p-cresol sink reaction (further metabolism unclear especially relevant for anaerobic conditions)
  mod <- sybil::addReact(mod, id="DM_cpd01042_c0", reactName="Sink needed for p-cresol", met="cpd01042[c0]", Scoef=-1, lb=0, ub=sybil::SYBIL_SETTINGS("MAXIMUM"), metComp = 1)
  mod@react_attr[which(mod@react_id == "DM_cpd01042_c0"),c("gs.origin","seed")] <- data.frame(gs.origin = 7, seed = "DM_cpd01042_c0", stringsAsFactors = F)

  
  mod <- add_missing_exchanges(mod) 
  
  # add metabolite & reaction attributes
  mod <- addMetAttr(mod, seed_x_mets = seed_x_mets)
  mod <- addReactAttr(mod)
  
  # add metabolite compartment list
  n.comp <- max(mod@met_comp, na.rm = T)
  if(n.comp == 2)
    mod@mod_compart <- c("c0","e0")
  if(n.comp == 3)
    mod@mod_compart <- c("c0","e0","p0")
  
  return(list(mod=mod, cand.rxns=dt.cand, rxn_x_genes=dt_genes))
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
source(paste0(script.dir,"/prepare_candidate_reaction_tables.R"))
source(paste0(script.dir,"/get_gene_logic_string.R"))
source(paste0(script.dir,"/addMetAttr.R"))
source(paste0(script.dir,"/addReactAttr.R"))
source(paste0(script.dir,"/parse_BMjson.R"))

# get options first
spec <- matrix(c(
  'blast.res', 'r', 1, "character", "Blast-results table generated by gapseq.sh.",
  'transporter.res', 't', 1, "character", "Blast-results table generated by transporter.sh.",
  'biomass', 'b', 2, "character", "Biomass reaction for network model. Options: Gram \"pos\" OR \"neg\" OR \"archaea\" OR \"auto\"? Default: \"auto\". Please note: if set to \"auto\", the external programms barrnap, blastn, and bedtools are required as well as argument \"-c\" for genome sequence.",
  'model.name', 'n', 2, "character", "Name of draft model network. Default: the basename of \"blast.res\"",
  'genome.seq', 'c', 2, "character", "If biomass \"-b\" is set to \"auto\", the genome sequence is required to search for 16S genes, which are used to predict the prokaryotic biomass reaction.",
  'high.evi.rxn.BS', "u", 2, "numeric", "Reactions with an associated blast-hit with a bitscore above this value will be added to the draft model as core reactions (i.e. high-sequence-evidence reactions)",
  'min.bs.for.core', "l", 2, "numeric", "Reactions with an associated blast-hit with a bitscore below this value will be considered just as reactions that have no blast hit.",
  'output.dir', 'o', 2, "character", "Directory to store results. Default: \".\" (alternatives not yet implemented)",
  'sbml.output', 's', 2, "logical", "Should the gapfilled model be saved as sbml? Default: FALSE (export not yet implemented)",
  'pathway.pred', 'p', 2, "character", "Pathway-results table generated by gapseq.sh.",
  'curve.alpha', 'a', 2, "numeric", "Exponent coefficient for transformation of bitscores to reaction weights for gapfilling. (Default: 1 (neg-linear))"
), ncol = 5, byrow = T)

opt <- getopt(spec)

# Help Screen
if ( is.null(opt$blast.res) | (is.null(opt$biomass) & is.null(opt$genome.seq)) | (!is.null(opt$biomass) && opt$biomass=="auto" && is.null(opt$genome.seq))) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# Setting defaults if required
if ( is.null(opt$model.name) ) { opt$model.name = NA_character_ }
if ( is.null(opt$output.dir) ) { opt$output.dir = "." }
if ( is.null(opt$sbml.output) ) { opt$sbml.output = F }
if ( is.null(opt$high.evi.rxn.BS) ) { opt$high.evi.rxn.BS = 200 }
if ( is.null(opt$min.bs.for.core) ) { opt$min.bs.for.core = 50 }
if ( is.null(opt$biomass) ) { opt$biomass = "auto" }
if ( is.null(opt$pathway.pred) ) { opt$pathway.pred = NA }
if ( is.null(opt$curve.alpha) ) { opt$curve.alpha = 1 }

# Arguments:
blast.res         <- opt$blast.res
transporter.res   <- opt$transporter.res
model.name        <- opt$model.name
biomass           <- opt$biomass
output.dir        <- opt$output.dir
genome.seq        <- opt$genome.seq
high.evi.rxn.BS   <- opt$high.evi.rxn.BS
min.bs.for.core   <- opt$min.bs.for.core
pathway.pred      <- opt$pathway.pred
curve.alpha       <- opt$curve.alpha

if(is.na(model.name)){
  model.name <- gsub("-all-Reactions.tbl","",basename(blast.res), fixed = T)
}

# construct draft model
mod <- build_draft_model_from_blast_results(blast.res = blast.res, 
                                            transporter.res = transporter.res,
                                            biomass = biomass, 
                                            model.name = model.name, 
                                            genome.seq = genome.seq, 
                                            high.evi.rxn.BS = high.evi.rxn.BS,
                                            min.bs.for.core = min.bs.for.core,
                                            script.dir = script.dir,
                                            pathway.pred = pathway.pred,
                                            curve.alpha = curve.alpha)

# save draft model and reaction weights and rxn-gene-table
saveRDS(mod$mod,file = paste0(model.name, "-draft.RDS"))
saveRDS(mod$cand.rxns,file = paste0(model.name, "-rxnWeights.RDS"))
saveRDS(mod$rxn_x_genes,file = paste0(model.name, "-rxnXgenes.RDS"))
# Write SBML
source(paste0(script.dir,"/sbml_write.R"))
write_gapseq_sbml(mod$mod, paste0(model.name,"-draft"))
