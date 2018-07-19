# This scripts uses three input files
# . dat/seed_reactions.tsv (original seed reaction database)
# . dat/corrections_seed_reactionDB.tsv (database of corrections to the original seed database)
# . dat/importantCIreactions.lst (list of charge-unbalanced reactions that are anyway accepted for the (full)-model)
correct_seed_rxnDB <- function(script.path) {
  require(data.table)
  require(stringr)
  require(sybil)
  source(paste0(script.dir, "/generate_rxn_stoich_hash.R"))
  
  mseed <- fread(paste0(script.dir,"/../dat/seed_reactions.tsv"), header=T, stringsAsFactors = F, 
                 colClasses = c("character","character","character","character","character","numeric",
                                "character","character","character","character","numeric","character","character",
                                "character","numeric","numeric",
                                "character","character","numeric","character","character"), na.strings = c("NA","null"))
  mseed[, gs.hash := generate_rxn_stoich_hash(stoichiometry, direction)]
  mseed.obs <- copy(mseed[is_obsolete == 1])
  mseed <- mseed[(status %in% c("OK","CPDFORMERROR") & is_obsolete == 0) | grepl("^CI:",status) | grepl("^OK.*\\|CI",status) | grepl("^MI:",status)]
  
  # duplicates with different reversibility as linked reactions (via hash)
  obs.map <- copy(mseed.obs[,.(id,is_obsolete,gs.hash,linked_reaction,name,abbreviation)])
  obs.map <- merge(obs.map,mseed[,.(id,gs.hash)],by="gs.hash", all.x = T)
  obs.map <- obs.map[is.na(id.y) & gs.hash != "ERR"]
  get.linked.rxn <- function(x) {
    cands <- unlist(strsplit(x,split = ";"))
    cands <- cands[cands %in% mseed$id]
    if(length(cands)>0)
      return(cands[1])
    return(NA)
  }
  obs.map[,id.y := get.linked.rxn(linked_reaction),by=id.x]
  obs.map <- merge(obs.map[,.(id.y,id.x,name,abbreviation,linked_reaction,is_obsolete)],
                   mseed[,-c("name","abbreviation","linked_reaction","is_obsolete")],
                   by.x="id.y",by.y="id")
  colnames(obs.map)[1:2] <- c("is_copy_of","id")
  dup.cor <- copy(obs.map)
  
  # duplicates with identical reversibility in mseed table
  obs.map <- copy(mseed.obs[,.(id,is_obsolete,gs.hash,linked_reaction)])
  obs.map <- merge(obs.map,mseed[,.(id,gs.hash)],by="gs.hash")
  obs.map <- obs.map[!duplicated(id.x)]
  obs.map <- obs.map[,.(id.x,is_obsolete,id.y)]
  colnames(obs.map) <- c("id","is_obsolete","is_copy_of")
  mseed.obs <- merge(mseed.obs,obs.map[,.(id,is_copy_of)],by="id")
  mseed[,is_copy_of := NA]
  mseed <- mseed[!(id %in% dup.cor$id) & !(id %in% mseed.obs$id)]
  mseed <- rbind(mseed,mseed.obs)
  mseed <- rbind(mseed,dup.cor)
  mseed[,gapseq.status:= "not.assessed"]
  mseed <- mseed[order(id)]
  
  # apply correction
  mseed.corr <- fread(paste0(script.dir,"/../dat/corrections_seed_reactionDB.tsv"), header= T, stringsAsFactors = F)
  if(system(paste0("grep \\\"\\\" ",script.path,"/../dat/corrections_seed_reactionDB.tsv"), ignore.stdout = T)==0) {
    stop("Error in seed-reactions corrections file: Double quotes!")
  }
  mseed.corr[is.na(rm.rxn), rm.rxn := FALSE]
  for(i in 1:nrow(mseed.corr)) {
    # Change is reversibility
    if(mseed.corr[i, new.rev]!="" & !is.na(mseed.corr[i, new.rev])) {
      mseed[id==mseed.corr[i, rnx.id] | is_copy_of==mseed.corr[i, rnx.id], 
            reversibility := mseed.corr[i, new.rev]]
      tmp <- mseed[id==mseed.corr[i, rnx.id], stoichiometry] # stoichiometry of corrected reaction
      mseed[is_copy_of==mseed.corr[i, rnx.id], 
            stoichiometry := tmp]
      mseed[id==mseed.corr[i, rnx.id] | is_copy_of==mseed.corr[i, rnx.id], 
            gapseq.status := "corrected"]
    }
    # change in stoichiometry
    if(mseed.corr[i, new.stoichiometry]!="" & !is.na(mseed.corr[i, new.stoichiometry])) {
      mseed[id==mseed.corr[i, rnx.id] | is_copy_of==mseed.corr[i, rnx.id],
            stoichiometry := mseed.corr[i, new.stoichiometry]]
      tmp <- mseed[id==mseed.corr[i, rnx.id], reversibility] # reversibility of corrected reaction
      mseed[is_copy_of==mseed.corr[i, rnx.id], 
            reversibility := tmp]
      mseed[id==mseed.corr[i, rnx.id] | is_copy_of==mseed.corr[i, rnx.id], 
            gapseq.status := "corrected"]
    }
    # remove reaction?
    if(mseed.corr[i, rm.rxn]==TRUE) {
      mseed[id==mseed.corr[i, rnx.id] | is_copy_of==mseed.corr[i, rnx.id], 
            gapseq.status := "removed"]
    }
    # New reaction?
    if(!(mseed.corr[i, rnx.id] %in% mseed$id) & mseed.corr[i, rnx.id] != "" & !is.na(mseed.corr[i, rnx.id])) {
      dt.n <- data.table(id            = mseed.corr[i, rnx.id],
                         abbreviation  = mseed.corr[i, abbreviation],
                         name          = mseed.corr[i, name],
                         stoichiometry = mseed.corr[i, new.stoichiometry],
                         is_transport  = 0,
                         reversibility = mseed.corr[i, new.rev],
                         direction     = mseed.corr[i, new.rev],
                         deltag        = 10000000,
                         deltagerr     = 10000000,
                         status        = "OK",
                         is_obsolete   = 0,
                         notes         = "Gapseq new",
                         gapseq.status = "corrected")
      mseed <- rbind(mseed,dt.n, fill=T)
    }
  }
  mseed[reversibility=="?", reversibility := direction]
  
  n.ind <- which(mseed$id=="rxn16318") # This is the reaction ID until which all reactions were tested/corrected to ensure that theres no futile cycle, that
                                       # would produce ATP from ADP + Pi without the consumtion of any energy source.
  tmp <- mseed[1:n.ind, id]
  mseed[id %in% tmp & gapseq.status == "not.assessed", gapseq.status := "approved"]
  mseed[is_obsolete==1 & is_copy_of %in% tmp & gapseq.status == "not.assessed", gapseq.status := "approved"]
  # importants reactions which are unfortunately charge unbalanced. Charge-unbalanced reactions are usually removed from database. However,
  # some are essential for metabolism and occure in Kbase-draft reconstructions.
  ci.rxns <- unlist(strsplit(readLines(paste0(script.dir,"/../dat/importantCIreactions.list"),warn = F)," "))
  mseed[(grepl("^CI:",status) | grepl("^OK.*\\|CI",status) | grepl("^MI:",status)) & gapseq.status == "approved" & !(id %in% ci.rxns),
        gapseq.status := "not.assessed"]
  table(mseed$gapseq.status)
  
  # re-generate definition and equation string for corrected reactions
  mseed[gapseq.status=="corrected",equation := getRxnEquaFromStoich(stoichiometry, reversibility, str.type = "equation")]
  mseed[gapseq.status=="corrected",definition := getRxnEquaFromStoich(stoichiometry, reversibility, str.type = "definition")]
  
  # Export of corrected database
  mseed <- mseed[,-"gs.hash"]
  fwrite(mseed, file=paste0(script.path,"/../dat/seed_reactions_corrected.tsv"), sep = "\t", quote = F)
}

getRxnEquaFromStoich <- function(stoichiometry, reversibility, str.type="definition") {
  out.strings <- c()
  for(i in 1:length(stoichiometry)) {
    rxn.info <- str_split(unlist(str_split(string = stoichiometry[i],pattern = ";")), pattern = ":", simplify = T)
    
    met.comp  <- rxn.info[,3]
    
    # stoichiometry
    met.scoef <- as.numeric(rxn.info[,1])
    
    # Ids for "equation"
    met.ids   <- paste0(rxn.info[,2],"[",met.comp,"]")
    met.ids.o <- paste0("(",abs(met.scoef),") ",met.ids)
    
    # Met names for reaction "definition"
    met.name  <- str_replace_all(rxn.info[,5],"\\\"","")
    met.name  <- paste0(met.name,"[",met.comp,"]")
    ind <- which(met.name=="" | is.na(met.name))
    met.name[ind] <- met.ids[ind]
    met.name.o <- paste0("(",abs(met.scoef),") ",met.name)
    
    # Reaction reversibility
    rev.sign <- ifelse(reversibility[i] == "<","<=",ifelse(reversibility[i] == ">","=>","<=>"))
    
    # LHS
    ind        <- which(met.scoef < 0)
    l.id.str   <- paste0(met.ids.o[ind], collapse = " + ")
    l.name.str <- paste0(met.name.o[ind], collapse = " + ")
    
    # RHS
    ind        <- which(met.scoef > 0)
    r.id.str   <- paste0(met.ids.o[ind], collapse = " + ")
    r.name.str <- paste0(met.name.o[ind], collapse = " + ")
    
    # final equation/definition assembly
    id.str   <- paste(l.id.str,rev.sign,r.id.str)
    name.str <- paste(l.name.str,rev.sign,r.name.str)
    
    if(str.type == "definition")
      out.strings <- c(out.strings, name.str)
    if(str.type == "equation")
      out.strings <- c(out.strings, id.str)
  }
  return(out.strings)
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

correct_seed_rxnDB(script.path = script.dir)