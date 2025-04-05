library(getopt)
suppressMessages(library(stringr))
suppressMessages(library(stringi))
suppressMessages(library(R.utils))
library(methods)
suppressMessages(library(IRanges))
options(error=traceback)

# test data
#blast.res <- "~/test/transl/HRGMv2_0240-all-Reactions.tbl"
#gram <- "pos"
#biomass <- "pos"
#genome.seq <- NA
#high.evi.rxn.BS <- 200
#transporter.res <- "~/test/MAG5-files/MAG5_protein-Transporter.tbl"

build_draft_model_from_blast_results <- function(blast.res, transporter.res, biomass = "auto", model.name = NA, genome.seq = NA,
                                                 script.dir, high.evi.rxn.BS = 200, pathway.pred = NA, min.bs.for.core = 50) {
  suppressMessages(require(data.table))
  setDTthreads(1)
  suppressMessages(require(stringr))
  suppressMessages(require(cobrar))

  # if biomass is "auto" get the taxonomy/Gram info from the input files
  findheader <- readLines(blast.res, n = 3)
  taxonomy <- str_match(findheader[3], "taxonomy=([^;]+)")[,2]
  gram <- str_match(findheader[3], "gram=([^;]+)")[,2]
  if(biomass == "auto") {
    if(length(taxonomy) != 1 || taxonomy %notin% c("Bacteria","neg","pos","Archaea")) {
      warning("Biomass could not be infered automatically from the input files. Assuming Bacteria.")
      biomass <- "Bacteria"
    } else {
      biomass <- taxonomy
    }
  }
  # If Bacteria, it still needs to be decided if Gram positive or negative:
  if(biomass == "Bacteria") {
    if(length(gram) == 1 && gram %in% c("neg","pos")) {
      biomass <- gram
    } else {
      source(paste0(script.dir,"/gram_by_network.R"))
      cat("Trying to predict biomass by metabolic network similarity\n")
      biomass <- find_gram_by_network(pathway.pred, script.dir)
    }
  }

  if(is.na(model.name))
    model.name <- gsub("-[a-z]+-Reactions.tbl","",basename(blast.res))

  if(biomass %in% c("auto","Bacteria","bacteria")) {
    forced_bacterial <- biomass %in% c("Bacteria", "bacteria")

    if(genome_mode == "nucl") {
      biomass <- system(paste0(script.dir,"/./predict_biomass_from16S.sh ",genome.seq), intern = T)

      # If user specified Bacteria but sequence-based prediction tool says archaeal
      if(biomass == "Archaea" & forced_bacterial)
        biomass <- "ambiguous"

      cat("\nPredicted biomass: ",biomass,"\n")
      if(biomass == "ambiguous" & !is.na(pathway.pred)) {
        cat("Trying to predict biomass by metabolic network similarity\n")
        biomass <- find_gram_by_network(pathway.pred, script.dir)
        cat("New predicted biomass: ",biomass,"\n")
      }

      if(!biomass %in% c("pos","neg", "archaea", "Archaea", "Gram_pos", "Gram_neg")) {
        stop("ERROR: Gram-staining prediction failed or ambiguous result. Please check whether genome sequence contains 16S rRNA gene(s).")

      }
    } else if(genome_mode == "prot") {
      # first predict domain if needed
      if(biomass %in% c("Auto","auto")) {
        hmmres <- tempfile()
        hmmprofiles <- R.utils::gunzip(paste0(script.dir,"/../dat/seq/hmm/domain.hmm.gz"),
                                       temporary = TRUE, overwrite = TRUE, remove = FALSE)
        system(paste0("hmmsearch --tblout ",hmmres," ",hmmprofiles," ",
                      genome.seq," > /dev/null"))
        biomass <- system(paste0("Rscript ",script.dir,"/predict_domain.R ",script.dir," ", hmmres), intern = TRUE)
        cat("Predicted domain:", biomass, "\n")
        tmp <- file.remove(c(hmmres, hmmprofiles))
        rm(tmp)
      }

      # second predict gram staining
      if(biomass %in% c("Bacteria", "bacteria")) {
        hmmres <- tempfile()
        hmmprofiles <- R.utils::gunzip(paste0(script.dir,"/../dat/seq/hmm/gram.hmm.gz"),
                                       temporary = TRUE, overwrite = TRUE, remove = FALSE)
        system(paste0("hmmsearch --tblout ",hmmres," ",hmmprofiles," ",
                      genome.seq," > /dev/null"))
        biomass <- system(paste0("Rscript ",script.dir,"/predict_gramstaining.R ",script.dir," ", hmmres), intern = TRUE)
        cat("Predicted Gram-staining:", biomass, "\n")
      }

    } else {
      biomass <- find_gram_by_network(pathway.pred, script.dir)
    }

  }

  seed_x_ec <- fread(paste0(script.dir,"/../dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv"), header=T)
  seed_x_name <- fread(paste0(script.dir,"/../dat/seed_reactions_corrected.tsv"), header=T, stringsAsFactors = F)
  #seed_x_metCyc <- fread(paste0(script.dir,"/../dat/mnxref_seed-other.tsv"), header = T)
  seed_x_aliases <- fread(paste0(script.dir,"/../dat/seed_Enzyme_Name_Reactions_Aliases.tsv"), header=T)
  seed_x_mets   <- fread(paste0(script.dir,"/../dat/seed_metabolites_edited.tsv"), header=T, stringsAsFactors = F, na.strings = c("null","","NA"))

  # Prepare candidate reaction tables for draft network and gapfilling
  dt.cand.tmp <- prepare_candidate_reaction_tables(blast.res, transporter.res, high.evi.rxn.BS, min.bs.for.core)
  dt      <- dt.cand.tmp$dt
  dt.cand <- dt.cand.tmp$dt.cand

  dt[!is.na(complex) & is.na(complex.status) & !pathway.status %in% c("full","threshold","keyenzyme"), complex.status := 0] # incomplete complexes

  mseed <- seed_x_name
  mseed <- mseed[gapseq.status %in% c("approved","corrected")]
  mseed <- mseed[order(id)]

  # Get all reactions / transporters that have either high sequence evidence (bitscore)
  dt_seed_single_and_there <- copy(dt[(bitscore >= high.evi.rxn.BS & status != "bad_blast") | pathway.status %in% c("full","threshold","keyenzyme")])
  dt_seed_single_and_there <- dt_seed_single_and_there[complex.status != 0 | is.na(complex.status)]
  dt_seed_single_and_there <- dt_seed_single_and_there[order(seed,-bitscore)]
  dt_seed_single_and_there <- dt_seed_single_and_there[!duplicated(seed)]

  # check if sequence names match conventions and are unique. if not, assign new ones.
  contig.names.full  <- unique(c(dt$stitle,dt.cand$stitle))
  contig.names.full  <- contig.names.full[contig.names.full != ""]
  n.contigs          <- length(contig.names.full)

  # Sequence title "corrections"
  # (1) cut everything after the first "space" occurrence, incl. the space
  # (2) Replace special characters
  # (3) Truncate sequence ID length if necessary
  # (4) Check if IDs are still unique, if not assign new unique sequence names
  # (5) Change sequence titles for gene names in model object

  new.stitle <- gsub(" .*$","", contig.names.full) # (1)
  new.stitle <- gsub("\\(|\\)|\\+|,|/|\\:|\\||\\&|\\[|\\]|\\{|\\}|'|\"","_", new.stitle) # (2)
  # (3)
  ind.toolong <- which(str_length(new.stitle) > 40)
  if(length(ind.toolong) > 0) {
    warning(paste0(length(ind.toolong)," genome contig/gene/protein sequence titles are very long. E.g.:\n\"",
                   contig.names.full[ind.toolong[1]], "\"\n",
                   "These titles are truncated to a maximum length of 40 characters. ",
                   "Please consider to adjust sequence titles to shorter identifiers before running gapseq.\n"))
    for(itl in ind.toolong)
      new.stitle[itl] <- substr(new.stitle[itl],1,40)
  }
  # (4)
  if(length(unique(new.stitle)) != n.contigs) {
    warning("Sequence identifiers do not match NCBI standard:\nhttps://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp\nAssigning new sequence names (\"sequence_[i]\").")
    new.stitle <- paste0("sequence_",1:n.contigs)
  }
  # (5)
  dt.new.ids <- data.table(old.contig = contig.names.full, new.contig = new.stitle) # dictionary
  dt.new.ids <- dt.new.ids[!duplicated(old.contig)]
  setkey(dt.new.ids, "old.contig")
  dt[stitle != "" & !is.na(stitle), stitle := dt.new.ids[stitle,new.contig]]
  dt.cand[stitle != "" & !is.na(stitle), stitle := dt.new.ids[stitle,new.contig]]

  # create gene list and attribute table
  cat("Creating Gene-Reaction list... ")

  dt_genes <- copy(dt[!is.na(bitscore)])
  #saveRDS(dt_genes, "~/dt_genes.RDS")
  dt_genes[, gene := gsub(" .*","",stitle)]
  dt_genes <- dt_genes[!duplicated(paste0(stitle, gene, seed, complex, sep = "$"))]

  dt_genes[grepl("Subunit undefined", complex), complex := "ZSubunit undefined"]
  dt_genes <- dt_genes[order(stitle,seed,complex,-bitscore)]
  dt_genes <- dt_genes[!duplicated(paste(stitle,seed,gene,sep = "$"))]
  dt_genes[grepl("ZSubunit undefined", complex), complex := "Subunit undefined"]

  dt_genes[, rm := FALSE]
  cat(length(unique(dt_genes[rm == F & seed %in% mseed$id, paste(gene, sep="$")])),
      "unique genes\n")

  # create subsys list and attribute table
  dt_subsys <- copy(dt[!is.na(pathway)])
  dt_subsys[, pathway := gsub("^\\||\\|$","",pathway)]
  subsys_unique <- unique(dt_subsys$pathway)

  rxns <- unique(dt_seed_single_and_there[,seed])

  # table of reactions for draft network
  if(any(mseed$id %in% rxns)) {
    mseed <- mseed[(id %in% rxns)]
  } else {
    mseed <- mseed[id %in% c("rxn00062")] # Adding ATP phosphohydrolase as dummy if no other reaction has been found. What kind of minimalist monster is this?
    warning("No reactions have been found to be added to the draft network. Returning a minimalist draft network.")
  }

  #remove duplicate reactions
  mseed[, rxn.hash := generate_rxn_stoich_hash(stoichiometry, reversibility)]
  mseed <- mseed[order(id, rxn.hash)]
  mseed <- mseed[!duplicated(rxn.hash)]

  cat("Constructing draft model: \n")
  mod <- new("ModelOrg", mod_name = model.name, mod_id = model.name)
  mod <- addCompartment(mod, c("c0","e0","p0"),
                        c("Cytosol","Extracellular space","Periplasm"))

  # add gapseq version info to model object
  gapseq_version <- system(paste0(script.dir,"/.././gapseq -v"), intern = T)[1]
  blast.header <- str_match(readLines(blast.res, n=2),"# Sequence DB md5sum: .*")
  if( any(!is.na(blast.header)) ){
      mod@mod_desc <- paste0(gapseq_version,"; ", na.omit(gsub("# ","",blast.header)))
  }else mod@mod_desc <- gapseq_version

  mod@react_attr$gs.origin <- integer(0L)
  mod@react_attr$name <- character(0L)
  mod@react_attr$ec <- character(0L)
  mod@react_attr$tc <- character(0L)
  mod@react_attr$exception <- integer(0L)
  mod@react_attr$complex.status <- integer(0L)

  reactAttrAdHocCols <- intersect(colnames(mod@react_attr),
                                  colnames(dt_seed_single_and_there))

  mod <- addSubsystem(mod, subsys_unique)
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

    ind.notnew.mets <- which(met.ids %in% mod@met_id)
    if(length(ind.notnew.mets) > 0)
      met.name[ind.notnew.mets] <- NA_character_

    # get reaction-associated stretches of DNA
    dtg.tmp <- dt_genes[seed == mseed[i,id] & bitscore >= high.evi.rxn.BS & rm == F, .(complex,gene)]
    setorder(dtg.tmp, gene) # ensure gene name order
    #if("Subunit undefined" %in% dtg.tmp$complex & !all(dtg.tmp$complex == "Subunit undefined")){
    #  dtg.tmp <- dtg.tmp[complex != "Subunit undefined"] # do not consider undefined subunits for gpr when defined subunits were found
    #}
    if(nrow(dtg.tmp)>0) {
      gpr.tmp <- get_gene_logic_string(dtg.tmp$complex, dtg.tmp$gene)
    } else {
      gpr.tmp <- ""
    }

    # get reaction-associated subsystems
    dts.tmp <- dt_subsys[seed == mseed[i,id], pathway]
    dts.tmp <- unique(dts.tmp)

    mod <- addReact(model = mod,
                    id = paste0(mseed[i,id],"_c0"),
                    met = met.ids,
                    Scoef = met.scoef,
                    metComp = mod@mod_compart[as.integer(met.comp)+1],
                    ub = ifelse(only.backwards, 0, COBRAR_SETTINGS("MAXIMUM")),
                    lb = ifelse(is.rev, -COBRAR_SETTINGS("MAXIMUM"), 0),
                    reactName = mseed[i, name],
                    metName = met.name,
                    gprAssoc = gpr.tmp,
                    subsystem = dts.tmp)
    if(mseed[i,id] %in% dt_seed_single_and_there[,seed]) {
      mod@react_attr[which(mod@react_id == paste0(mseed[i,id],"_c0")),reactAttrAdHocCols] <- as.data.frame(dt_seed_single_and_there[seed == mseed[i,id]])[1,reactAttrAdHocCols]
    }

  }

  mod@react_attr$gs.origin <- 0
  mod@react_attr$gs.origin[mod@react_attr$bitscore < high.evi.rxn.BS] <- 9 # Added due to Pathway Topology criteria
  #mod <- add_reaction_from_db(mod, react = c("rxn13782","rxn13783","rxn13784"), gs.origin = 6) # Adding pseudo-reactions for Protein biosynthesis, DNA replication and RNA transcription
  mod <- add_missing_diffusion(mod)

  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
  # Adding transporters if pertinent pathway/reaction is present #
  # GS origin code: 5                                            #
  #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

  # in case of electron bifurcating butanoyl-CoA dehydrogenase (NAD+, ferredoxin) presence:
  # add butyrate transporter as well
  if(any(grepl("rxn90001", mod@react_id)) & all(!grepl("rxn05683", mod@react_id))) {
    mod <- add_reaction_from_db(mod, react = "rxn05683", gs.origin = 5)
  }

  # in case of tryptophan degradation to IPA (indole-3-propionate) add IPA transport
  if(any(grepl("rxn43343", mod@react_id)) & any(grepl("rxn45361", mod@react_id)) & any(grepl("rxn00483", mod@react_id)) & any(grepl("rxn01447", mod@react_id)) & all(!grepl("rxn90116", mod@react_id))) {
    mod <- add_reaction_from_db(mod, react = "rxn90116", gs.origin = 5)
  }

  # in case of phenylalanine degradation to PPA (phenylpropanoate) add PPA transport
  if(any(grepl("rxn00493", mod@react_id)) & any(grepl("rxn00997", mod@react_id)) & any(grepl("rxn07603", mod@react_id)) & any(grepl("rxn40746", mod@react_id)) & all(!grepl("rxn09182", mod@react_id))) {
    mod <- add_reaction_from_db(mod, react = "rxn09182", gs.origin = 5)
  }

  # in case of Tyrosine degradation to Phloretate add Phloretate transport
  if(any(grepl("rxn00527", mod@react_id)) & any(grepl("rxn02393", mod@react_id)) & any(grepl("rxn46948", mod@react_id)) & any(grepl("rxn46031", mod@react_id)) & all(!grepl("rxn90117", mod@react_id))) {
    mod <- add_reaction_from_db(mod, react = "rxn90117", gs.origin = 5)
  }

  cat("\n")

  # Adding Biomass reaction
  if(biomass == "neg" | biomass == "Gram_neg"){
    ls.bm <- parse_BMjson(paste0(script.dir, "/../dat/biomass/biomass_Gram_neg.json"), seed_x_mets)
    dt.bm <- ls.bm$bmS
  } else if(biomass == "pos" | biomass == "Gram_pos"){
    ls.bm <- parse_BMjson(paste0(script.dir, "/../dat/biomass/biomass_Gram_pos.json"), seed_x_mets)
    dt.bm <- ls.bm$bmS
  } else if(biomass == "archaea" | biomass == "Archaea"){
    ls.bm <- parse_BMjson(paste0(script.dir, "/../dat/biomass/biomass_archaea.json"), seed_x_mets)
    dt.bm <- ls.bm$bmS
  } else if(file.exists(biomass)){
    cat("Using user-defined biomass definition...\n")
    ls.bm <- parse_BMjson(biomass, seed_x_mets)
    dt.bm <- ls.bm$bmS
  } else {
    stop("Invalid value for biomass option '-b'.")
  }
  # mod@mod_attr <- data.frame(annotation = paste0("tax_domain:",ls.bm$domain))
  mod@mod_attr$annotation <- paste0("tax_domain:",ls.bm$domain)

  if(biomass %in% c("neg","pos","Gram_neg","Gram_pos")) {

    # remove menaquinone8 from anaerobic Biomass if de novo biosynthesis pathway is absent

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

  mnames <- paste0(dt.bm$name,"-",dt.bm$comp,"0")
  mnames <- ifelse(dt.bm$id %in% mod@met_id, NA_character_, dt.bm$name)

  mod <- addReact(mod,id = "bio1",
                  met = dt.bm$id,
                  Scoef = dt.bm$Scoef,
                  lb = 0, ub = COBRAR_SETTINGS("MAXIMUM"),
                  obj = 1,
                  reactName = ls.bm$name,
                  metName = mnames,
                  metComp = mod@mod_compart[ifelse(dt.bm$comp=="c",1,ifelse(dt.bm$comp=="e",2,3))],
                  SBOTerm = "SBO:0000629")
  mod@react_attr[which(mod@react_id == "bio1"),c("gs.origin","seed")] <- data.frame(gs.origin = 6, seed = "bio1", stringsAsFactors = F)

  # add p-cresol sink reaction (further metabolism unclear especially relevant for anaerobic conditions)
  mod <- addReact(mod, id="DM_cpd01042_c0", reactName="Sink needed for p-cresol",
                  met="cpd01042[c0]", metName="p-Cresol-c0", Scoef=-1, lb=0,
                  ub=COBRAR_SETTINGS("MAXIMUM"), metComp = 1)
  mod@react_attr[which(mod@react_id == "DM_cpd01042_c0"),c("gs.origin","seed")] <- data.frame(gs.origin = 7, seed = "DM_cpd01042_c0", stringsAsFactors = F)


  mod <- add_missing_exchanges(mod)

  # add metabolite & reaction attributes
  mod <- addMetAttr(mod, seed_x_mets = seed_x_mets)
  mod <- addReactAttr(mod)
  mod <- addGeneAttr(mod, dt_genes)

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
source(paste0(script.dir,"/addGeneAttr.R"))
source(paste0(script.dir,"/parse_BMjson.R"))

# get options first
spec <- matrix(c(
  'blast.res', 'r', 1, "character", "Blast-results table generated by gapseq find.",
  'transporter.res', 't', 1, "character", "Blast-results table generated by gapseq find-transport.",
  'biomass', 'b', 2, "character", "Biomass reaction for network model. Options: Gram \"pos\" OR \"neg\" OR \"archaea\" OR \"auto\" OR path to user-defined biomass reaction. Default: \"auto\". See Details.",
  'model.name', 'n', 2, "character", "Name of draft model network. Default: the basename of \"blast.res\"",
  'genome.seq', 'c', 2, "character", "If biomass \"-b\" is set to \"auto\", the genome sequence is required to predict the prokaryotic biomass reaction.",
  'high.evi.rxn.BS', "u", 2, "numeric", "Reactions with an associated blast-hit with a bitscore above this value will be added to the draft model as core reactions (i.e. high-sequence-evidence reactions)",
  'min.bs.for.core', "l", 2, "numeric", "Reactions with an associated blast-hit with a bitscore below this value will be considered just as reactions that have no blast hit.",
  'output.dir', 'f', 2, "character", "Path to directory, where output files will be saved (default: current directory)",
  'depr.output.dir', 'o', 2, "character", "deprecated. Use flag\"-f\" instead",
  'sbml.no.output', 's', 2, "logical", "Do not save draft model as sbml file. Default: Save as SBML",
  'pathway.pred', 'p', 2, "character", "Pathway-results table generated by gapseq find.",
  'help' , 'h', 0, "logical", "help"
), ncol = 5, byrow = T)

opt <- getopt(spec)

# Help Screen
if ( !is.null(opt$help) | is.null(opt$blast.res) | (is.null(opt$biomass) & is.null(opt$genome.seq)) | (!is.null(opt$biomass) && opt$biomass=="auto" && is.null(opt$genome.seq))) {

  if(opt$biomass=="auto" && is.null(opt$genome.seq)) {
    cat("\nInput error: No genome supplied for biomass (bacteria / archaea) prediction. Please provide a genome sequence with the option \"-c|--genome.seq\" or specify a biomass catagory (e.g. \"Archaea\",\"Gram_neg\",\"Gram_pos\") with option\"-b|--biomass\".\n\n")
  }

  cat(getopt(spec, usage=TRUE))

  cat("\n")
  cat("Details:\n")
  cat("\"biomass\" (-b): See gapseq documentation (https://gapseq.readthedocs.io/en/latest/database/biomassReaction.html) for more details.\n")
  q(status=1)
}

# Setting defaults if required
if ( is.null(opt$model.name) ) { opt$model.name = NA_character_ }
if ( is.null(opt$genome.seq) ) { opt$genome.seq = NA_character_ }
if ( is.null(opt$output.dir) ) { opt$output.dir = "." }
if ( is.null(opt$sbml.no.output) ) { opt$sbml.no.output = F } else { opt$sbml.no.output = T }
if ( is.null(opt$high.evi.rxn.BS) ) { opt$high.evi.rxn.BS = 200 }
if ( is.null(opt$min.bs.for.core) ) { opt$min.bs.for.core = 50 }
if ( is.null(opt$biomass) ) { opt$biomass = "auto" }
if ( is.null(opt$pathway.pred) ) { opt$pathway.pred = NA }

# deprecation notice for flag '-o'
if(!is.null(opt$depr.output.dir)) {
  warning("Deprecation notice: Flag '-o' is now replaced with flag '-f'. Please adjust your script(s), as the flag will be removed in a future release.")
  if(is.null(opt$output.dir))
    opt$output.dir <- opt$depr.output.dir
}

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
sbml.no.output    <- opt$sbml.no.output

# create output directory if not already there
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

if(is.na(model.name)){
  model.name <- gsub("-[a-z]+-Reactions\\.tbl$|-[a-z]+-Reactions\\.tbl\\.gz$","",
                     basename(blast.res),
                     fixed = FALSE)
}

if(!(biomass %in% c("auto","Auto","pos","neg","archaea","Archaea","bacteria",
                   "Bacteria", "Gram_pos", "Gram_neg") || file.exists(biomass))) {
  stop("Invalid value for biomass option '-b'.")
} else if (!dir.exists(output.dir) || file.access(output.dir, mode = 2) == -1)  {
  stop(paste("Output directory",output.dir,"cannot be created or is not writable."))
} else {
  # construct draft model
  mod <- build_draft_model_from_blast_results(blast.res = blast.res,
                                              transporter.res = transporter.res,
                                              biomass = biomass,
                                              model.name = model.name,
                                              genome.seq = genome.seq,
                                              high.evi.rxn.BS = high.evi.rxn.BS,
                                              min.bs.for.core = min.bs.for.core,
                                              script.dir = script.dir,
                                              pathway.pred = pathway.pred)

  # remove empty subsystems
  subsRm <- which(apply(mod$mod@subSys,2,FUN = function(x) all(x==FALSE)))
  mod$mod <- cobrar::rmSubsystem(mod$mod, subsRm)

  # save draft model and reaction weights and rxn-gene-table
  saveRDS(mod$mod,file = paste0(output.dir,"/",model.name, "-draft.RDS"))
  saveRDS(mod$cand.rxns,file = paste0(output.dir,"/",model.name, "-rxnWeights.RDS"))
  saveRDS(mod$rxn_x_genes,file = paste0(output.dir,"/",model.name, "-rxnXgenes.RDS"))

  # Write SBML
  if(!sbml.no.output) {
    source(paste0(script.dir,"/sbml_write.R"))
    write_gapseq_sbml(mod$mod, paste0(output.dir,"/",model.name,"-draft"))
  }
}
