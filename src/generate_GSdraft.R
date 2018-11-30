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
# topology.evidence <- "../CarveMe-tests/faa/GCF_000005845.2_ASM584v2_genomic-Transporter.lst"

# Please note: If Game-Staining is set to "auto" it requires the external programms barrnap, usearch and bedtools
build_draft_model_from_blast_results <- function(blast.res, topo.evi = NA, gram = "auto", model.name = NA, genome.seq = NA, script.dir, high.evi.rxn.BS = 200) {
  suppressMessages(require(data.table))
  suppressMessages(require(stringr))
  #suppressMessages(require(sybil))
  
  if(is.na(model.name))
    model.name <- gsub("-all-Reactions.tbl","",basename(blast.res), fixed = T)
  
  if(is.na(topo.evi[1]))
    topo.evi <- c()
  
  if(gram=="auto") {
    if(grepl("\\.gz$", genome.seq)) {
      suppressMessages(require(R.utils))
      genome.seq <- gunzip(genome.seq, remove = F, temporary = T, overwrite = T)
    } else {
      
    }
    system(paste0("cat ",genome.seq, " | sed '/^[[:space:]]*$/d' | tr -d '\15\32' > ", genome.seq, ".tmp")) # remove empty lines => problems with bedtools
    system(paste0("barrnap --quiet ",genome.seq, ".tmp | grep 16S > ", genome.seq, ".gff"))
    system(paste0("bedtools getfasta -s -name+ -fi ",genome.seq,".tmp -bed ", genome.seq,".gff -fo ", genome.seq, ".16S.fasta"))
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
  
  
  dt <- fread(blast.res, header=T)
  dt <- dt[,.(rxn, name, ec, qseqid, pident, evalue, bitscore, qcovs, stitle, sstart, send, pathway, status, pathway.status, dbhit)]
  #dt[str_count(ec,"\\.")==1, ec := paste0(ec,".-.-")] # TODO: RM
  #dt[str_count(ec,"\\.")==2, ec := paste0(ec,".-")] # TODO: RM
  
  # A specific fix for the issue with reactions 1.3.8.1 and 1.3.8.13
  calc_seq_overlap <- function(astart, aend, bstart, bend) {
    return(length(intersect(astart:aend,bstart:bend))/((length(astart:aend)+length(bstart:bend))/2))
  }
  
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
     dt <- dt[-rm.ids]
  }

  
  # # Add those reactions that are going to be added due to pathway completion
  # if(!is.na(pwy.pred)) {
  #   dt.pwy <- fread(pwy.pred, header=T, stringsAsFactors = F)
  #   pres.pwys <- dt.pwy[Prediction==TRUE, ID] # all present pathways
  #   
  #   # next, get all reactions associated with present pathways
  #   dt.pwy1 <- fread(paste0(script.dir,"/../dat/meta_pwy.tbl"), header=T, stringsAsFactors = F)
  #   dt.pwy2 <- fread(paste0(script.dir,"/../dat/custom_pwy.tbl"), header=T, stringsAsFactors = F)
  #   dt.pwy <- rbind(dt.pwy1, dt.pwy2, fill = T)
  #   dt.pwy <- dt.pwy[id %in% pres.pwys]
  #   
  #   dt.pwy.topo.rea <- data.table(rxn = character(0),name = character(0), ec = character(0), topo.evidence = character(0))
  #   for(i in 1:nrow(dt.pwy)) {
  #     tmp.rxn  <- unlist(str_split(dt.pwy[i, reaId],","))
  #     name.str <- dt.pwy[i, reaName]
  #     name.str <- str_replace_all(name.str, "&Delta;","&Delta$")
  #     tmp.name <- unlist(str_split(name.str,";"))
  #     tmp.name <- str_replace_all(tmp.name,"&Delta$", "&Delta;")
  #     tmp.ec   <- unlist(str_split(dt.pwy[i, reaEc],","))
  #     
  #     if(length(tmp.rxn) != length(tmp.name) | length(tmp.rxn) != length(tmp.ec))
  #       stop(paste("Error: Number ob EC numbers, reactions IDs and/or reaction names not identical for Pathway:",dt.pwy[i,id]))
  #     
  #     dt.tmp <- data.table(rxn = tmp.rxn, name = tmp.name, ec = tmp.ec, topo.evidence = paste("Topology Evidence:",dt.pwy[i,id]))
  #     dt.pwy.topo.rea <- rbind(dt.pwy.topo.rea, dt.tmp)
  #   }
  #   
  #   dt.pwy.topo.rea <- dt.pwy.topo.rea[!duplicated(paste(rxn,name,sep="$$"))]
  #   
  #   dt <- merge(dt, dt.pwy.topo.rea[,.(rxn,topo.evidence)], by = "rxn", all.x = T)
  #   
  #   dt <- rbind(dt, dt.pwy.topo.rea[!(rxn %in% dt[,rxn])], fill=T)
  # } else {
  #   dt[,topo.evidence := NA_character_]
  # }
  # 
  # 
  # # 1. assign modelSEED ID via ec-number
  # dt_1 <- merge(dt, seed_x_ec[,.(`External ID`,seed=`MS ID`)], by.x = "ec", by.y = "External ID")
  # 
  # # 2. assign modelSEED ID via mnxref namespace
  # dt_2 <- merge(dt, seed_x_metCyc[,.(seed, other)] , by.x = "rxn", by.y = "other")
  # 
  # # for remaining unmapped MetaCyc reaction try this:
  # # 3. assign modelSEED ID via reaction name
  # matches <- unique(c(dt_1[,paste(rxn, name, sep = "$")], dt_2[,paste(rxn, name, sep = "$")]))
  # dt_3 <- merge(dt[!(paste(rxn,name,sep="$") %in% matches)], seed_x_name[,.(seed=id,name)] , by.x = "name", by.y = "name")
  # matches <- c(matches, dt_3[,paste(rxn, name, sep = "$")])
  # dt_4 <- merge(dt[!(paste(rxn,name,sep="$") %in% matches)], seed_x_aliases[,.(seed=seed.rxn.id,enzyme.name)], by.x = "name",by.y = "enzyme.name")
  # matches <- c(matches, dt_4[,paste(rxn, name, sep = "$")])
  # 
  # dt_seed <- rbindlist(list(dt_1,dt_2,dt_3,dt_4))
  # dt_seed <- splitcol2rows_mget(dt_seed, "seed", "|")
  # 
  # #TODO: for developing purposes (will be removed later on.) Write file with unmapped meta-cyc reactions (to seed DB)
  # dt_unmap <- dt[!(paste(rxn,name,sep="$") %in% matches), .(rxn, name,ec, bitscore)]
  # fwrite(dt_unmap, file=paste0(blast.res,"_unmapps_MCrxns.csv"), quote = FALSE, sep="\t")
  
  names(dt)[names(dt) == "dbhit"] = "seed" 
  dt <- splitcol2rows_mget(dt, "seed", " ")
  
  mseed <- seed_x_name
  mseed <- mseed[gapseq.status %in% c("approved","corrected")]
  mseed <- mseed[order(id)]
  
  # Get all reactions that have either high sequence evidence (bitscore) or pathway-topology evidence (pathway status)
  dt_seed_single_and_there <- copy(dt[bitscore >= high.evi.rxn.BS | pathway.status %in% c("full","treshold","keyenzyme")])
  dt_seed_single_and_there <- dt_seed_single_and_there[order(seed,-bitscore)]
  dt_seed_single_and_there <- dt_seed_single_and_there[!duplicated(seed)]
  
  rxns <- unique(dt_seed_single_and_there[,seed])
  
  mseed <- mseed[(id %in% rxns) | (id %in% topo.evi)]
  #mseed <- mseed[(id %in% rxns)]
  
  #remove duplicate reactions
  mseed[, rxn.hash := generate_rxn_stoich_hash(stoichiometry, reversibility)]
  mseed <- mseed[order(id, rxn.hash)]
  mseed <- mseed[!duplicated(rxn.hash)]
  
  mod <- sybil::modelorg(name = model.name, id = model.name)
  mod@react_attr <- data.frame(rxn = character(0), name = character(0), ec = character(0), qseqid = character(0), pident = numeric(0), evalue = numeric(0),
                               bitscore = numeric(0), qcovs = numeric(0), stitle = character(0), sstart = numeric(0), send = numeric(0),
                               pathway = character(0), status = character(0), pathway.status = character(0), seed = character(0), stringsAsFactors = F)
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
    
    mod <- sybil::addReact(model = mod, 
                    id = paste0(mseed[i,id],"_c0"), 
                    met = met.ids,
                    Scoef = met.scoef,
                    reversible = is.rev, 
                    metComp = as.integer(met.comp)+1,
                    ub = ifelse(only.backwards, 0, 1000),
                    lb = ifelse(is.rev, -1000, 0),
                    reactName = mseed[i, name], 
                    metName = met.name)
    if(mseed[i,id] %in% dt_seed_single_and_there[,seed])
      mod@react_attr[which(mod@react_id == paste0(mseed[i,id],"_c0")),] <- as.data.frame(dt_seed_single_and_there[seed == mseed[i,id]])
  }
  mod@react_attr$gs.origin <- 0
  mod@react_attr$gs.origin[mod@react_attr$bitscore < high.evi.rxn.BS] <- 9 # Added due to Pathway Topology criteria
  mod@react_attr$gs.origin[is.na(mod@react_attr$bitscore) & is.na(mod@react_attr$topo.evidence)] <- 5 # Transporter
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
  'topology.evidence', 't', 2, "character", "File with space-separated list of reactions (SEED-namespace), that have gapseq-topology evidence.",
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
if ( is.null(opt$topology.evidence) ) { opt$topology.evidence = NA_character_ }

# Arguments:
blast.res         <- opt$blast.res
model.name        <- opt$model.name
gram              <- opt$gram
output.dir        <- opt$output.dir
genome.seq        <- opt$genome.seq
high.evi.rxn.BS   <- opt$high.evi.rxn.BS
topology.evidence <- opt$topology.evidence

if(is.na(model.name))
  model.name <- gsub("-all-Reactions.tbl","",basename(blast.res), fixed = T)

# loading transporters and reactions with topology evidence
topo.evi <- c()
for( file in unlist(str_split(topology.evidence, ",")) ){
  input.rxns <- readLines(file)
  topo.evi <- c(topo.evi, unlist(str_split(input.rxns, " ")))
}
topo.evi <- unique(topo.evi)
cat("\nLoading input reactions:", length(topo.evi), "\n")

# construct draft model
mod <- build_draft_model_from_blast_results(blast.res = blast.res, 
                                            topo.evi = topo.evi,
                                            gram = gram, 
                                            model.name = model.name, 
                                            genome.seq = genome.seq, 
                                            high.evi.rxn.BS = high.evi.rxn.BS,
                                            script.dir = script.dir)

# save draft model and reaction weights
saveRDS(mod$mod,file = paste0(model.name, ".RDS"))
saveRDS(mod$cand.rxns,file = paste0(model.name, "-rxnWeights.RDS"))

