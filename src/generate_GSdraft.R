library(getopt)
suppressMessages(library(stringr))
suppressMessages(library(stringi))
library(methods)
options(error=traceback)

# test data
# blast.res <- "../CarveMe-tests/faa/GCF_000005845.2_ASM584v2_genomic-all-Reactions.tbl"
# gram <- "auto"
# genome.seq <- "../CarveMe-tests/faa/GCF_000005845.2_ASM584v2_genomic.fna.gz"
# high.evi.rxn.BS <- 200
# transporter.res <- "../CarveMe-tests/faa/GCF_000005845.2_ASM584v2_genomic-Transporter.tbl"

# Please note: If Game-Staining is set to "auto" it requires the external programms barrnap, usearch and bedtools
build_draft_model_from_blast_results <- function(blast.res, transporter.res, gram = "auto", model.name = NA, genome.seq = NA, script.dir, high.evi.rxn.BS = 200) {
  suppressMessages(require(data.table))
  suppressMessages(require(stringr))
  #suppressMessages(require(sybil))
  
  if(is.na(model.name))
    model.name <- gsub("-all-Reactions.tbl","",basename(blast.res), fixed = T)
  
  if(gram=="auto") {
    if(grepl("\\.gz$", genome.seq)) {
      suppressMessages(require(R.utils))
      genome.seq <- gunzip(genome.seq, remove = F, temporary = T, overwrite = T)
    } else {
      
    }
    system(paste0("cat ",genome.seq, " | sed '/^[[:space:]]*$/d' | tr -d '\15\32' > ", genome.seq, ".tmp")) # remove empty lines => problems with bedtools
    system(paste0("barrnap --quiet ",genome.seq, ".tmp | grep 16S > ", genome.seq, ".gff"))
    system(paste0("bedtools getfasta -s -name -fi ",genome.seq,".tmp -bed ", genome.seq,".gff -fo ", genome.seq, ".16S.fasta"))
    #system(paste0("rnammer -S bac -m ssu -f ",genome.seq,".16S.fasta ",genome.seq))
    #system(paste0("usearch -sinaps ",genome.seq,".16S.fasta -db " ,script.dir,"/../dat/seq/Bacteria/16S_graminfo/16S_gramposinfo.fna -attr grampos -tabbedout ",genome.seq,".graminfo.tsv -strand plus"))
    
    gram <- system(paste0(script.dir,"/./predict_Gram_staining_from16S.sh ",genome.seq,".16S.fasta"), intern = T)
    
    #gram.dt <- fread(paste0(genome.seq,".graminfo.tsv"), header = F)
    #gram.dt[, max.score := max(V3)]
    #gram.dt <- gram.dt[V3 > max.score * .95]
    #gram <- names(sort(table(gram.dt$V2),decreasing=TRUE)[1])
    cat("\nPredicted gram staining: ",gram,"\n")
    if(gram == "ambiguous") {
      stop("ERROR: Gram-staining prediction failed or ambiguous result. Please check whether genome sequence contains 16S rRNA gene(s).")
    }
    
    # clean up
    system(paste0("rm ", genome.seq, ".tmp"))
    system(paste0("rm ", genome.seq, ".tmp.fai"))
    system(paste0("rm ", genome.seq, ".gff"))
    system(paste0("rm ", genome.seq, ".16S.fasta"))
  }
  
  seed_x_ec <- fread(paste0(script.dir,"/../dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv"), header=T)
  seed_x_name <- fread(paste0(script.dir,"/../dat/seed_reactions_corrected.tsv"), header=T, stringsAsFactors = F)
  seed_x_metCyc <- fread(paste0(script.dir,"/../dat/mnxref_seed-other.tsv"), header = T)
  seed_x_aliases <- fread(paste0(script.dir,"/../dat/seed_Enzyme_Name_Reactions_Aliases.tsv"), header=T)
  
  # Read reaction blast results
  dt <- fread(blast.res, header=T, stringsAsFactors = F)
  dt <- dt[,.(rxn, name, ec, tc = NA_character_, qseqid, pident, evalue, bitscore, qcovs, stitle, sstart, send, pathway, status, pathway.status, seed = dbhit)]
  # Read transporter blast results
  dt.trans <- fread(transporter.res, header=T, stringsAsFactors = F)
  dt.trans <- dt.trans[,.(rxn = id, name = paste("transport",tc,sub,sep="-"), ec = NA_character_, tc, qseqid, pident, evalue, bitscore, qcovs, stitle, 
                          sstart, send, pathway = NA_character_, status = NA_character_, pathway.status = NA_character_, seed = rea)]
  dt.trans[bitscore >= high.evi.rxn.BS, status := "good_blast"]
  dt.trans[bitscore <  high.evi.rxn.BS, status := "bad_blast"]
  
  # a little helper function to calculate the overlap of two genetic regions
  calc_seq_overlap <- function(astart, aend, bstart, bend) {
    out <- numeric(length = length(bstart))
    for(i in seq_along(bstart)){
      out[i] <- length(intersect(astart:aend,bstart[i]:bend[i]))/((length(astart:aend)+length(bstart[i]:bend[i]))/2)
    }
    return(out)
  }
  
  # A specific fix for the issue with reactions 1.3.8.1 and 1.3.8.13
  if("1.3.8.1" %in% dt$ec & "1.3.8.13" %in% dt$ec) {
     rm.ids <- c()
     one.hits.id <- which(dt$ec == "1.3.8.1" & !is.na(dt$bitscore))
     two.hits.id <- which(dt$ec == "1.3.8.13" & !is.na(dt$bitscore))
     if(length(one.hits.id)>0 & length(two.hits.id)>0) {
       for(i in one.hits.id) {
         for(j in two.hits.id) {
           ol <- calc_seq_overlap(dt[i,sstart],dt[i,send],dt[j,sstart],dt[j,send])
           if(ol > 0.2) {
             if(dt[i, bitscore] > max(dt[two.hits.id, bitscore]))
               rm.ids <- c(rm.ids, j)
             if(dt[j, bitscore] > max(dt[one.hits.id, bitscore]))
               rm.ids <- c(rm.ids, i)
           }
         }
       }
     }
     rm.ids <- unique(rm.ids)
     if(length(rm.ids > 0))
       dt <- dt[-rm.ids]
  }

  # prepare reaction/transporter blast results table
  dt <- splitcol2rows_mget(dt, "seed", " ")
  dt.trans <- splitcol2rows_mget(dt.trans, "seed", ",")
  dt <- rbind(dt, dt.trans)
  
  mseed <- seed_x_name
  mseed <- mseed[gapseq.status %in% c("approved","corrected")]
  mseed <- mseed[order(id)]
  
  # Get all reactions / transporters that have either high sequence evidence (bitscore)
  # TODO: Here we are currently loosing the information of isozymes
  dt_seed_single_and_there <- copy(dt[bitscore >= high.evi.rxn.BS | pathway.status %in% c("full","treshold","keyenzyme")])
  dt_seed_single_and_there <- dt_seed_single_and_there[order(seed,-bitscore)]
  dt_seed_single_and_there <- dt_seed_single_and_there[!duplicated(seed)]
  
  # create gene list and attribute table
  cat("Creating Gene-Reaction list... \n")
  dt_genes <- copy(dt[!is.na(bitscore)])
  dt_genes[, gene := paste(sstart,send, sep = ":")]
  dt_genes <- dt_genes[!duplicated(paste(stitle,seed,gene,sep = "$"))]
  dt_genes <- dt_genes[order(stitle,seed,-bitscore)]
  dt_genes[, rm := F]
  dt_genes$itmp <- 1:nrow(dt_genes)
  
  for(i in 1:nrow(dt_genes)) {
    astart <- dt_genes[i, sstart]
    aend <- dt_genes[i, send]
    astitle <- dt_genes[i, stitle]
    aseed <- dt_genes[i, seed]
    dt_g_tmp <- dt_genes[i < itmp & rm==F & stitle == astitle & seed == aseed]
    if(nrow(dt_g_tmp)>0) {
      dt_g_tmp[calc_seq_overlap(astart, aend, sstart, send) > 0.5, rm := T]
      ind.rm <- dt_g_tmp[rm==T, itmp]
      dt_genes[itmp %in% ind.rm, rm := T]
    }
  }
  
  
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
  mod@react_attr <- data.frame(rxn = character(0), name = character(0), ec = character(0), tc = character(0), qseqid = character(0),
                               pident = numeric(0), evalue = numeric(0), bitscore = numeric(0), qcovs = numeric(0),
                               stitle = character(0), sstart = numeric(0), send = numeric(0), pathway = character(0),
                               status = character(0), pathway.status = character(0), seed = character(0), stringsAsFactors = F)
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
    
    # get reaction-associated strechtes of dna 
    # TODO: handle mutli-contig genomes
    dtg.tmp <- dt_genes[seed == mseed[i,id] & bitscore >= high.evi.rxn.BS, gene]
    dtg.tmp <- unique(dtg.tmp)
    if(length(dtg.tmp > 0))
      gpr.tmp <- paste(dtg.tmp, collapse = " | ")
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
                    ub = ifelse(only.backwards, 0, 1000),
                    lb = ifelse(is.rev, -1000, 0),
                    reactName = mseed[i, name], 
                    metName = met.name,
                    gprAssoc = gpr.tmp,
                    subSys = dts.tmp)
    if(mseed[i,id] %in% dt_seed_single_and_there[,seed])
      mod@react_attr[which(mod@react_id == paste0(mseed[i,id],"_c0")),] <- as.data.frame(dt_seed_single_and_there[seed == mseed[i,id]])
  }
  mod@react_attr$gs.origin <- 0
  mod@react_attr$gs.origin[mod@react_attr$bitscore < high.evi.rxn.BS] <- 9 # Added due to Pathway Topology criteria
  mod <- add_reaction_from_db(mod, react = c("rxn13782","rxn13783","rxn13784"), gs.origin = 6) # Adding pseudo-reactions for Protein biosynthesis, DNA replication and RNA transcription
  #mod@genes_table <- copy(dt[bitscore>0])
  cat("\n")
  
  # Adding Biomass reaction
  if(gram == "neg")
    dt.bm <- fread(paste0(script.dir, "/../dat/seed_biomass.DT_gramNeg.tsv"))
  if(gram == "pos")
    dt.bm <- fread(paste0(script.dir, "/../dat/seed_biomass.DT_gramPos.tsv"))
  
  mod <- sybil::addReact(mod,id = "bio1", met = dt.bm$id, Scoef = dt.bm$stoich, reversible = F, lb = 0, ub = 1000, obj = 1, 
                  reactName = paste0("Biomass reaction ",ifelse(gram=="neg","(gram -)","(gram +)")), metName = dt.bm$name, metComp = dt.bm$comp)
  mod@react_attr[which(mod@react_id == "bio1"),c("gs.origin","seed")] <- data.frame(gs.origin = 6, seed = "bio1", stringsAsFactors = F)
  
  mod <- add_missing_exchanges(mod) 
  
  #
  # construct gapfill candidate reactions and their weight
  #
  dt.cand <- dt[bitscore < high.evi.rxn.BS]
  dt.cand[, max.bs := max(bitscore), by = "seed"]
  dt.cand <- dt.cand[max.bs == bitscore]
  dt.cand <- dt.cand[!duplicated(seed)]
  dt.cand[, max.bs := NULL]
  dt.cand[, weight := 1 - (bitscore / high.evi.rxn.BS)]
  
  return(list(mod=mod, cand.rxns=dt.cand))
}

splitcol2rows_mget <- function(dtInput, col2split, sep){
  dtInput <- dtInput[, .(tmp.add.col = unlist(strsplit(get(col2split),sep,T))), by=names(dtInput)]
  
  dtInput[, c(col2split):=NULL];
  setnames(dtInput, 'tmp.add.col', col2split); 
  return(dtInput);
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
  'blast.res', 'r', 1, "character", "Blast-results table generated by gapseq.sh.",
  'transporter.res', 't', 1, "character", "Blast-results table generated by transporter.sh.",
  'gram', 'g', 2, "character", "Gram \"pos\" OR \"neg\" OR \"auto\"? Default: \"auto\". Please note: if set to \"auto\", the external programms barrnap, usearch, and bedtools are required.",
  'model.name', 'n', 2, "character", "Name of draft model network. Default: the basename of \"blast.res\"",
  'genome.seq', 'c', 2, "character", "If gram is set to \"auto\", the genome sequence is required to search for 16S genes, which are used to predict gram-staining.",
  'high.evi.rxn.BS', "b", 2, "numeric", "Reactions with an associated blast-hit with a bitscore above this value will be added to the draft model as core reactions (i.e. high-sequence-evidence reactions)",
  'output.dir', 'o', 2, "character", "Directory to store results. Default: \".\" (alternatives not yet implemented)",
  'sbml.output', 's', 2, "logical", "Should the gapfilled model be saved as sbml? Default: FALSE (export not yet implemented)"
), ncol = 5, byrow = T)

opt <- getopt(spec)

# Help Screen
if ( is.null(opt$blast.res) | (is.null(opt$gram) & is.null(opt$genome.seq)) | (!is.null(opt$gram) && opt$gram=="auto" && is.null(opt$genome.seq))) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# Setting defaults if required
if ( is.null(opt$model.name) ) { opt$model.name = NA_character_ }
if ( is.null(opt$output.dir) ) { opt$output.dir = "." }
if ( is.null(opt$sbml.output) ) { opt$sbml.output = F }
if ( is.null(opt$high.evi.rxn.BS) ) { opt$high.evi.rxn.BS = 200 }
if ( is.null(opt$gram) ) { opt$gram = "auto" }

# Arguments:
blast.res         <- opt$blast.res
transporter.res   <- opt$transporter.res
model.name        <- opt$model.name
gram              <- opt$gram
output.dir        <- opt$output.dir
genome.seq        <- opt$genome.seq
high.evi.rxn.BS   <- opt$high.evi.rxn.BS

if(is.na(model.name))
  model.name <- gsub("-all-Reactions.tbl","",basename(blast.res), fixed = T)

# construct draft model
mod <- build_draft_model_from_blast_results(blast.res = blast.res, 
                                            transporter.res = transporter.res,
                                            gram = gram, 
                                            model.name = model.name, 
                                            genome.seq = genome.seq, 
                                            high.evi.rxn.BS = high.evi.rxn.BS,
                                            script.dir = script.dir)

# save draft model and reaction weights
saveRDS(mod$mod,file = paste0(model.name, ".RDS"))
saveRDS(mod$cand.rxns,file = paste0(model.name, "-rxnWeights.RDS"))
